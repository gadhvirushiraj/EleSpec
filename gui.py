import tkinter as tk
import customtkinter as ctk
import matplotlib.pyplot as plt
from typing import Union, Callable
from const import elnames, molnames
from detect import driver_gui
import os
from collections import OrderedDict
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class GUI():
    search_items = OrderedDict()
    numeric_search_items = OrderedDict()
    items_is_ele = None

    def __init__(self) -> None:

        # system settings
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # configure ctk
        self.root = ctk.CTk()
        self.root.geometry("900x700")
        self.root.minsize(900,700)
        self.root.title('EleSpec')

        self.add_plot(None, None)
        self.add_left_bar()
        self.add_right_frame()

        self.root.mainloop()

    def add_plot(self, fig, ax):

        # configure plot style
        if True:
            color = 'k'
            bg = 'white'
        else:
            color = 'white'
            bg = '#232323'

        # configure plot 
        if fig == None or ax == None :
            fig, ax = plt.subplots(figsize=(self.root.winfo_width()/100, self.root.winfo_height()/200))

        ax.set_facecolor('white') # change background color to dark grey -> #232323
        fig.patch.set_facecolor('white') # change background color to dark grey -> #232323

        ax.tick_params(axis='x', colors=color)
        ax.tick_params(axis='y', colors=color)
        plt.yticks(fontsize=5)
        plt.xticks(fontsize=5)

        ax.spines[['right', 'top']].set_visible(False)
        ax.spines['left'].set_color(color)
        ax.spines['bottom'].set_color(color)

        ax.margins(y=0.1)
        ax.ticklabel_format(useOffset=False, style='scientific', axis='y')
        plt.ylim(ax.get_ylim()[0]*0.5, ax.get_ylim()[1]*1.7)
        plt.xlabel('Wavelength (Å)', fontsize=5)
        plt.ylabel('Flux', fontsize=5)
    
        fig.subplots_adjust(bottom=0.2)

        # add plot to gui
        canvas = FigureCanvasTkAgg(fig, master=self.root)
        canvas.draw()
        canvas.get_tk_widget().place(relx=0, rely=0, relwidth=1, relheight=1/2)

        # add toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.root)
        toolbar.update()
        toolbar.place(relx=0, rely=1/2, relwidth=1, relheight=0.05)

    def add_left_bar(self):

        self.frame1 = ctk.CTkFrame(self.root)
        self.frame1.place(relx=0, rely=1/2+0.05, relwidth=1, relheight=1-1/2-0.05)

        # add file path label
        ctk.CTkLabel(self.frame1, text="File Path (absolute)", font=("Arial", 15)).place(relx=0.03, rely=0.02)

        # add text entry 
        self.text_box = ctk.CTkEntry(self.frame1,width=200, font=("Arial", 15), placeholder_text='/Users/ME/mydata')
        self.text_box.place(relx=0.03, rely=0.1, relwidth=2/7)

        # to change dropdown menu values
        def edit_combo():
            if self.is_ele.get() == 1:
                combo.configure(values=elnames[:-1])
                combo.set('all-atoms')
            else:
                combo.configure(values=molnames)
                combo.set('all-molecules')

        # add selection radiobutton
        self.is_ele = tk.IntVar()
        self.is_ele.set(1)
        ele_check = ctk.CTkRadioButton(self.frame1, radiobutton_height=20, radiobutton_width=20 , text="Search Atom", font=("Arial", 15), variable=self.is_ele, value=1, command=edit_combo)
        mol_check = ctk.CTkRadioButton(self.frame1, radiobutton_height=20, radiobutton_width=20, text="Search Molecule", font=("Arial", 15), variable=self.is_ele, value=0, command=edit_combo)
        ele_check.place(relx=0.03, rely=0.22)
        mol_check.place(relx=0.03, rely=0.33)

        # adding info label
        ctk.CTkLabel(self.frame1, text="Select Item:", font=("Arial", 15)).place(relx=0.03, rely=0.42)

        # dropdown menu
        combo = ctk.CTkComboBox(self.frame1,values= elnames[:-1], font=("Arial", 15))
        combo.place(relx=0.03, rely=0.50, relwidth=1.5/7)

        # ADD button beside dropdown menu
        ctk.CTkButton(self.frame1, text="ADD", font=("Arial", 15),command= lambda: self.add_item(combo.get())).place(relx=0.03+1.5/7, rely=0.50, relwidth=0.5/7)

        # add info label
        ctk.CTkLabel(self.frame1, text="spectrum lowerbound", font=("Arial", 15)).place(relx=0.03, rely=0.60)

        # add lowerbound spinbox
        self.xmin = FloatSpinbox(self.frame1, width=150, step_size=50)
        self.xmin.place(relx=0.03, rely=0.68, relwidth=2/7)
        self.xmin.set(5880)

        # add info label
        ctk.CTkLabel(self.frame1, text="spectrum upperbound", font=("Arial", 15)).place(relx=0.03, rely=0.78)

        # add upperbound spinbox
        self.xmax= FloatSpinbox(self.frame1, width=150, step_size=50)
        self.xmax.place(relx=0.03, rely=0.86, relwidth=2/7)
        self.xmax.set(5930)
    
    def add_right_frame(self):

        # secoundary frame (placed in right bottom of root)
        self.frame2 = ctk.CTkFrame(self.frame1, fg_color="grey28", corner_radius=0)
        self.frame2.place(relx=0.35, rely=0, relwidth=0.65, relheight=1)

        # items to search
        ctk.CTkLabel(self.frame2, text="Items to search", font=("Arial", 15)).place(relx=0.03, rely=0.02)
        self.search_label = ctk.CTkLabel(self.frame2, text=",".join(self.search_items), font=("Arial", 15), fg_color='grey18', anchor="w", corner_radius=8)
        self.search_label.place(relx=0.03, rely=0.1, relwidth=0.78)

        # add clear button to right of search label
        ctk.CTkButton(self.frame2, text="Clear", font=("Arial", 15), command=self.clear_list).place(relx=0.82, rely=0.1, relwidth=0.15)

        # add info label
        ctk.CTkLabel(self.frame2, text="Del Atom Parameter", font=("Arial", 15)).place(relx=0.03, rely=0.22)

        # text box -> delatom
        self.delatom = ctk.CTkEntry(self.frame2,width=200, font=("Arial", 15),justify='right')
        self.delatom.place(relx=0.28, rely=0.22, relwidth=2/7)
        self.delatom.insert(0, '10')

        # add info label
        ctk.CTkLabel(self.frame2, text="Del Mol Parameter", font=("Arial", 15)).place(relx=0.03, rely=0.32)
        
        # text box -> delmol
        self.delmol = ctk.CTkEntry(self.frame2,width=200, font=("Arial", 15),justify='right')
        self.delmol.place(relx=0.28, rely=0.32, relwidth=2/7)
        self.delmol.insert(0, '0')

        # add info labbel
        ctk.CTkLabel(self.frame2, text="Representation Mode", font=("Arial", 15)).place(relx=0.6, rely=0.22)

        # add representaion mode (numeric, text, ifacode) radio button, placed at same horizontal level
        self.represent_mode = tk.IntVar()
        self.represent_mode.set(1)
        ctk.CTkRadioButton(self.frame2, radiobutton_height=20, radiobutton_width=20, text="Text", font=("Arial", 15), variable=self.represent_mode, value=1, command= lambda: self.switch(1)).place(relx=0.6, rely=0.32)
        ctk.CTkRadioButton(self.frame2, radiobutton_height=20, radiobutton_width=20, text="Numeric", font=("Arial", 15), variable=self.represent_mode, value=0, command= lambda: self.switch(0)).place(relx=0.75, rely=0.32)

        # add info labbel
        ctk.CTkLabel(self.frame2, text="Normalize Data", font=("Arial", 15)).place(relx=0.6, rely=0.42)

        # add normalize radio button
        self.norm = tk.IntVar()
        self.norm.set(1)
        ctk.CTkRadioButton(self.frame2, radiobutton_height=20, radiobutton_width=20, text="ON", font=("Arial", 15), variable=self.norm, value=1).place(relx=0.6, rely=0.52)
        ctk.CTkRadioButton(self.frame2, radiobutton_height=20, radiobutton_width=20, text="OFF", font=("Arial", 15), variable=self.norm, value=0).place(relx=0.75, rely=0.52)

        # add button to bottom of frame
        ctk.CTkButton(self.frame2, text="Process", font=("Arial", 15), command=self.process).pack(side="bottom", pady=20)


    def add_item(self, value):

        if self.items_is_ele == None:
            self.items_is_ele = self.is_ele.get()
        
        if self.is_ele.get() == 1 and self.items_is_ele == False:
            tk.messagebox.showerror("Error", "Atom cannot be added to molecule list")
            return
        
        if self.is_ele.get() == 0 and self.items_is_ele == True:
            tk.messagebox.showerror("Error", "Molecule cannot be added to atom list")
            return

        if value.strip() in ["all-atoms", 'all-molecules']:
            self.search_items = OrderedDict()
            self.numeric_search_items = OrderedDict()

            self.search_items[value.strip()] = None
            self.numeric_search_items[value.strip()] = None

            self.switch(1)

        elif 'all-atoms' in self.search_items or 'all-molecules' in self.search_items:
            tk.messagebox.showerror("Error", "all already in list, clear the list to add other items")
        else:
            self.search_items[value.strip()] = None

            if self.is_ele.get() == 1:
                self.numeric_search_items[str(elnames.index(value.strip()) * 100)] = None
            else:
                self.numeric_search_items[str(molnames.index(value.strip()))] = None

            self.switch(self.represent_mode.get())

    def clear_list(self):
        self.search_items = OrderedDict()
        self.numeric_search_items = OrderedDict()
        self.items_is_ele = None
        self.search_label.configure(text="")

    def switch(self, value):
        if value == 1:
            self.search_label.configure(text= "  " + ", ".join(self.search_items))
        elif value == 0:
            self.search_label.configure(text= "  " + ", ".join(self.numeric_search_items))
               
    def process(self):

        # check if file path is valid
        if not os.path.isfile(self.text_box.get()):
            tk.messagebox.showerror("Error", "File path is invalid")
            return
        
        # check if search item list is empty
        if list(self.search_items.keys()) == []:
            tk.messagebox.showerror("Error", "add atleast one item to search")
            return  
        
        # lowerbound < upperbound check
        if self.xmin.get() >= self.xmax.get():
            tk.messagebox.showerror("Error", "Lowerbound must be strictly smaller than upperbound")
            return
        
        # check if delatom and delmol is not empty
        if self.delatom.get() == "" or self.delmol.get() == "":
            tk.messagebox.showerror("Error", "Delatom and Delmol is needed")
            return
        
        # check if delatom and delmol is numeric integer
        if self.delatom.get().isnumeric() == False and self.delmol.get().isnumeric() == False and int(self.delatom.get()) < 0 and int(self.delmol.get()) < 0:
            tk.messagebox.showerror("Error", "Delatom and Delmol must be numeric integer")
            return
        
        
        search = []
        if self.represent_mode.get() == 1:
            search = list(self.search_items.keys())
        elif self.represent_mode.get() == 0:
            search = list(self.numeric_search_items.keys())
        
        
        fig, ax = driver_gui(self.text_box.get(),'k',self.norm.get(),search,self.xmin.get(),self.xmax.get(),self.represent_mode.get(),not self.is_ele.get(),int(self.delatom.get()),int(self.delmol.get()))
        self.add_plot(fig, ax)


class FloatSpinbox(ctk.CTkFrame):
    def __init__(self, *args,
                 width: int = 100,
                 height: int = 32,
                 step_size: Union[int, float] = 1,
                 command: Callable = None,
                 **kwargs):
        super().__init__(*args, width=width, height=height, **kwargs)

        self.step_size = step_size
        self.command = command

        self.configure(fg_color=("gray78", "gray28"))  # set frame color

        self.grid_columnconfigure((0, 2), weight=0)  # buttons don't expand
        self.grid_columnconfigure(1, weight=1)  # entry expands

        self.subtract_button = ctk.CTkButton(self, text="-", width=height-6, height=height-6,
                                                       command=self.subtract_button_callback)
        self.subtract_button.grid(row=0, column=0, padx=(3, 0), pady=3)

        self.entry = ctk.CTkEntry(self, width=width-(2*height), height=height-6, border_width=0)
        self.entry.grid(row=0, column=1, columnspan=1, padx=3, pady=3, sticky="ew")

        self.add_button = ctk.CTkButton(self, text="+", width=height-6, height=height-6,
                                                  command=self.add_button_callback)
        self.add_button.grid(row=0, column=2, padx=(0, 3), pady=3)

        # default value
        self.entry.insert(0, "0.0")

    def add_button_callback(self):
        if self.command is not None:
            self.command()
        try:
            value = float(self.entry.get()) + self.step_size
            self.entry.delete(0, "end")
            self.entry.insert(0, value)
        except ValueError:
            return

    def subtract_button_callback(self):
        if self.command is not None:
            self.command()
        try:
            value = float(self.entry.get()) - self.step_size
            self.entry.delete(0, "end")
            self.entry.insert(0, value)
        except ValueError:
            return

    def get(self) -> Union[float, None]:
        try:
            return float(self.entry.get())
        except ValueError:
            return None

    def set(self, value: float):
        self.entry.delete(0, "end")
        self.entry.insert(0, str(float(value)))


# add all to elnames and molnames
elnames[0] = "all-atoms"
molnames[0] = "all-molecules"

gui = GUI()