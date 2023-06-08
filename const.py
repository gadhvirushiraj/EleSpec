# Defining all the elements
elnames=['??','H','He','Li','Be','B','C','N','O','F',
'Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc',
'Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',
'Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh',
'Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La',
'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',
'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl',
'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np',
'Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db',
'Sg','Bh','Hs','Mt','Ds','Rg','??']


# Defining all the ionization stages
ionstages=['I','II','III','IV','V','VI','VII','VIII','IX','X',
'XI','XII','XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX',
'XXI','XXII','XXIII','XXIV','XXV','XXVI','XXVII','XXVIII','XXIX','XXX']

# Defining all the molecules
molnames = ["" for x in range(327)]
molnames[0] = '??????????'
molnames[1] = 'H2O'
molnames[2] = '46TiO'
molnames[3] = '47TiO'
molnames[4] = '48TiO'
molnames[5] = '49TiO'
molnames[6] = '50TiO'
molnames[7] = 'H2'
molnames[8] = '12CH'
molnames[9] = '13CH'
molnames[10] = '14NH'
molnames[11] = '15NH'
molnames[12] = '16OH'
molnames[13] = '17OH'
molnames[14] = '18OH'
molnames[15] = '24MgH'
molnames[16] = '25MgH'
molnames[17] = '26MgH'
molnames[18] = '28SiH'
molnames[19] = '29SiH'
molnames[20] = '30SiH'
molnames[21] = '12C12C'
molnames[22] = '12C13C'
molnames[23] = '13C13C'
molnames[24] = '12C14N'
molnames[25] = '13C14N'
molnames[26] = '12C15N'
molnames[27] = '13C15N'
molnames[28] = '12C16O'
molnames[29] = '13C16O'
molnames[30] = '12C17O'
molnames[31] = '13C17O'
molnames[32] = '12C18O'
molnames[33] = '13C18O'
molnames[34] = '28Si16O'
molnames[35] = '29Si16O'
molnames[36] = '30Si16O'
molnames[37] = '28Si18O'
molnames[38] = '46TiO [J]'
molnames[39] = '47TiO [J]'
molnames[40] = '48TiO [J]'
molnames[41] = '49TiO [J]'
molnames[42] = '50TiO [J]'
molnames[43] = 'H2-16O[H]'
molnames[44] = 'H2-18O[H]'
molnames[45] = 'H2-17O[H]'
molnames[46] = 'H16OD [H]'
molnames[47] = '12C16O2'
molnames[48] = '13C16O2'
molnames[49] = '12C16O18O'
molnames[50] = '12C16O17O'
molnames[51] = '13C16O18O'
molnames[52] = '13C16O17O'
molnames[53] = '12C18O2'
molnames[54] = '12C18O17O'
molnames[55] = '16O3'
molnames[56] = '16O16O18O'
molnames[57] = '16O18O16O'
molnames[58] = '14N2-16O'
molnames[59] = '14N15N16O'
molnames[60] = '15N14N16O'
molnames[61] = '14N2-18O'
molnames[62] = '14N2-17O'
molnames[63] = '12C16O[H]'
molnames[64] = '13C16O[H]'
molnames[65] = '12C18O[H]'
molnames[66] = '12C17O[H]'
molnames[67] = '13C18O[H]'
molnames[68] = '12CH4'
molnames[69] = '13CH4'
molnames[70] = '12CH3D'
molnames[71] = '16O2'
molnames[72] = '16O18O'
molnames[73] = '16O17O'
molnames[74] = '14N16O'
molnames[75] = '15N16O'
molnames[76] = '14N18O'
molnames[77] = '32S16O2'
molnames[78] = '34S16O2'
molnames[79] = '14N16O2'
molnames[80] = '14NH3'
molnames[81] = '15NH3'
molnames[82] = 'H14N16O3'
molnames[83] = '16OH [H]'
molnames[84] = '18OH [H]'
molnames[85] = '16OD [H]'
molnames[86] = 'H19F'
molnames[87] = 'H35Cl'
molnames[88] = 'H37Cl'
molnames[89] = 'H79Br'
molnames[90] = 'H81Br'
molnames[91] = 'H127I'
molnames[92] = '35Cl16O'
molnames[93] = '37Cl16O'
molnames[94] = '16O12C32S'
molnames[95] = '16O12C34S'
molnames[96] = '16O13C32S'
molnames[97] = '18O12C32S'
molnames[98] = 'H2-12C16O'
molnames[99] = 'H2-13C16O'
molnames[100] = 'H2-12C18O'
molnames[101] = 'H16O35Cl'
molnames[102] = 'H16O37Cl'
molnames[103] = '14N2'
molnames[104] = 'H12C14N'
molnames[105] = 'H13C14N'
molnames[106] = 'H12C15N'
molnames[107] = '12CH3-35Cl'
molnames[108] = '12CH3-37Cl'
molnames[109] = 'H2-16O2'
molnames[110] = '12C2H2'
molnames[111] = '12C13CH2'
molnames[112] = '12C2H6'
molnames[113] = '31PH3'
molnames[114] = '12C16O19F2'
molnames[115] = '32S19F6'
molnames[116] = 'H2-32S'
molnames[117] = 'H2-34S'
molnames[118] = 'H2-33S'
molnames[119] = '12C16O[G]'
molnames[120] = '13C16O[G]'
molnames[121] = '12C18O[G]'
molnames[122] = '12C17O[G]'
molnames[123] = '13C18O[G]'
molnames[124] = '13C17O[G]'
molnames[125] = '14C16O[G]'
molnames[126] = '12C14N[J]'
molnames[127] = '13C14N[J]'
molnames[128] = '28SiO [N]'
molnames[129] = '29SiO [N]'
molnames[130] = '30SiO [N]'
molnames[131] = 'H3+'
molnames[132] = 'YO'
molnames[133] = '90ZrO'
molnames[134] = '91ZrO'
molnames[135] = '92ZrO'
molnames[136] = '93ZrO'
molnames[137] = '94ZrO'
molnames[138] = '95ZrO'
molnames[139] = '96ZrO'
molnames[140] = 'H2O[AMES]'
molnames[141] = '12CH4[GH]'
molnames[142] = '13CH4[GH]'
molnames[143] = '12CH3D[GH]'
molnames[144] = '12CH4[GE]'
molnames[145] = '13CH4[GE]'
molnames[146] = '12CH3D[GE]'
molnames[147] = '46TiO [S]'
molnames[148] = '47TiO [S]'
molnames[149] = '48TiO [S]'
molnames[150] = '49TiO [S]'
molnames[151] = '50TiO [S]'
molnames[152] = 'FeH'
molnames[153] = 'H16OD [A]'
molnames[154] = '12CH4/STDS'
molnames[155] = '13CH4/STDS'
molnames[156] = '12CD4/STDS'
molnames[157] = '13CD4/STDS'
molnames[158] = '12CH4/STDS'
molnames[159] = '13CH4/STDS'
molnames[160] = '12CD4/STDS'
molnames[161] = '13CD4/STDS'
molnames[162] = 'VO[AMES]'
molnames[163] = 'CrH[AMES]'
molnames[164] = 'H2O[SCAN]'
molnames[165] = '12CH4/STDn'
molnames[166] = '12CH4/STDn'
molnames[167] = '12CH4/STDn'
molnames[168] = '12CH4/STDn'
molnames[169] = '12CH4/STDn'
molnames[170] = '12CH4/STDn'
molnames[171] = '12CH4/STDn'
molnames[172] = '12CH4/STDn'
molnames[173] = '12CH4/STDn'
molnames[174] = '12CH4/STDn'
molnames[175] = '12CH4/STDn'
molnames[176] = '12CH4/STDn'
molnames[177] = '12CH4/STDn'
molnames[178] = '12CH4/STDn'
molnames[179] = '12CH4/STDn'
molnames[180] = '12CH4/STDn'
molnames[181] = '12CH4/STDn'
molnames[182] = '13CH4/STDn'
molnames[183] = '13CH4/STDn'
molnames[184] = '13CH4/STDn'
molnames[185] = '13CH4/STDn'
molnames[186] = '13CH4/STDn'
molnames[187] = '24MgH[P+P]'
molnames[188] = '25MgH[P+P]'
molnames[189] = '26MgH[P+P]'
molnames[190] = '40CaH[P+P]'
molnames[191] = '44CaH[P+P]'
molnames[192] = 'H2O[BT]'
molnames[193] = '46TiO[Plez]'
molnames[194] = '47TiO[Plez]'
molnames[195] = '48TiO[Plez]'
molnames[196] = '49TiO[Plez]'
molnames[197] = '50TiO[Plez]'
molnames[198] = '51VO[Plez]'
molnames[199] = 'H2-16O'
molnames[200] = 'H2-18O'
molnames[201] = 'H2-17O'
molnames[202] = 'H16OD'
molnames[203] = 'H18OD'
molnames[204] = 'H17OD'
molnames[205] = '12C16O2'
molnames[206] = '13C16O2'
molnames[207] = '12C16O18O'
molnames[208] = '12C16O17O'
molnames[209] = '13C16O18O'
molnames[210] = '13C16O17O'
molnames[211] = '12C18O2'
molnames[211] = '12C18O2'
molnames[212] = '12C18O17O'
molnames[213] = '16O3'
molnames[214] = '16O16O18O'
molnames[217] = '16O17O16O'
molnames[215] = '16O18O16O'
molnames[216] = '16O16O17O'
molnames[218] = '14N2-16O'
molnames[219] = '14N15N16O'
molnames[220] = '15N14N16O'
molnames[221] = '14N2-18O'
molnames[222] = '14N2-17O'
molnames[223] = '12C16O [H]'
molnames[224] = '13C16O [H]'
molnames[225] = '12C18O [H]'
molnames[226] = '12C17O [H]'
molnames[227] = '13C18O [H]'
molnames[228] = '13C17O [H]'
molnames[229] = '12CH4'
molnames[230] = '13CH4'
molnames[231] = '12CH3D'
molnames[232] = '16O2'
molnames[233] = '16O18O'
molnames[234] = '16O17O'
molnames[235] = '14N16O'
molnames[236] = '15N16O'
molnames[237] = '14N18O'
molnames[238] = '32S16O2'
molnames[239] = '34S16O2'
molnames[240] = '14N16O2'
molnames[241] = '14NH3'
molnames[242] = '15NH3'
molnames[243] = 'H14N16O3'
molnames[244] = '16OH [H]'
molnames[245] = '18OH [H]'
molnames[246] = '16OD [H]'
molnames[247] = 'H19F'
molnames[248] = 'H35Cl'
molnames[249] = 'H37Cl'
molnames[250] = 'H79Br'
molnames[251] = 'H81Br'
molnames[252] = 'H127I'
molnames[253] = '35Cl16O'
molnames[254] = '37Cl16O'
molnames[255] = '16O12C32S'
molnames[256] = '16O12C34S'
molnames[257] = '16O13C32S'
molnames[258] = '16O12C33S'
molnames[259] = '18O12C32S'
molnames[260] = 'H2-12C16O'
molnames[261] = 'H2-13C16O'
molnames[262] = 'H2-12C18O'
molnames[263] = 'H16O35Cl'
molnames[264] = 'H16O37Cl'
molnames[265] = '14N2'
molnames[266] = 'H12C14N'
molnames[267] = 'H13C14N'
molnames[268] = 'H12C15N'
molnames[269] = '12CH3-35Cl'
molnames[270] = '12CH3-37Cl'
molnames[271] = 'H2-16O2'
molnames[272] = '12C2H2'
molnames[273] = '12C13CH2'
molnames[274] = '12C2H6'
molnames[275] = '31PH3'
molnames[276] = '12C16O19F2'
molnames[277] = '32S19F6'
molnames[278] = 'H2-32S'
molnames[279] = 'H2-34S'
molnames[280] = 'H2-33S'
molnames[281] = 'H12C16OOH'
molnames[282] = 'H16O2'
molnames[283] = '35ClONO2'
molnames[284] = '37ClONO2'
molnames[285] = '14N16O+'
molnames[286] = 'H16O79Br'
molnames[287] = 'H16O81Br'
molnames[288] = '12C12CH4'
molnames[289] = '12C13CH4'
molnames[290] = '12CH316OH'
molnames[291] = '40CaH[Plez'
molnames[292] = 'FeH [B]'
molnames[293] = 'TiH [B]'
molnames[294] = 'CrH [B]'
molnames[295] = '14NH3 [S] '
molnames[296] = '12CH4 [Br]'
molnames[297] = '12CH4 [Br]'
molnames[298] = '12CH4 [Br]'
molnames[299] = '12C16O2'
molnames[300] = '13C16O2'
molnames[301] = '12C16O18O'
molnames[302] = '12C16O17O'
molnames[303] = '13C16O18O'
molnames[304] = '13C16O17O'
molnames[305] = '12C18O2'
molnames[306] = '12C18O17O'
molnames[307] = '13C18O2'
molnames[308] = '16OH [H08]'
molnames[309] = '18OH [H08]'
molnames[310] = '16OD [H08]'
molnames[311] = '31PH3[H08]'
molnames[312] = '12C16O2'
molnames[313] = '13C16O2'
molnames[314] = '12C16O18O'
molnames[315] = '12C16O17O'
molnames[316] = '13C16O18O'
molnames[317] = '13C16O17O'
molnames[318] = '12C18O2'
molnames[319] = '12C16O2'
molnames[320] = '13C16O2'
molnames[321] = '12C16O18O'
molnames[322] = '12C16O17O'
molnames[323] = '13C16O18O'
molnames[324] = '13C16O17O'
molnames[325] = '12C18O2'
molnames[326] = '12C18O17O'


# translation table to France's codes

ifacod = [0 for x in range(327)]

#-- H2O:
ifacod[1] = 112 
#
#-- TiO:
ifacod[2] = 111 
ifacod[3] = 111 
ifacod[4] = 111 
ifacod[5] = 111 
ifacod[6] = 111 
#
#-- TiO [Jorgensen]:
ifacod[38] = 111 
ifacod[39] = 111 
ifacod[40] = 111 
ifacod[41] = 111 
ifacod[42] = 111 
#
#-- H2:
ifacod[7] = 113 
#
#-- CH:
ifacod[8] = 103 
ifacod[9] = 103 
#
#-- NH:
ifacod[10] = 104 
ifacod[11] = 104 
#
#-- OH:
ifacod[12] = 102 
ifacod[13] = 102 
ifacod[14] = 102 
#
#-- MgH:
ifacod[15] = 108 
ifacod[16] = 108 
ifacod[17] = 108 
#
#-- SiH:
ifacod[18] = 110 
ifacod[19] = 110 
ifacod[20] = 110 
#
#-- C2:
ifacod[21] = 105 
ifacod[22] = 105 
ifacod[23] = 105 
#
#-- CN:
ifacod[24] = 106 
ifacod[25] = 106 
ifacod[26] = 106 
ifacod[27] = 106 
#
#-- CO:
ifacod[28] = 107 
ifacod[29] = 107 
ifacod[30] = 107 
ifacod[31] = 107 
ifacod[32] = 107 
ifacod[33] = 107 
#
#-- SiO:
ifacod[34] = 121 
ifacod[35] = 121 
ifacod[36] = 121 
ifacod[37] = 121 
#
#-- H2O HITRAN :
ifacod[43] = 112 
ifacod[44] = 112 
ifacod[45] = 112 
ifacod[46] = 112 
#
#-- CO2 HITRAN :
ifacod[47] = 116 
ifacod[48] = 116 
ifacod[49] = 116 
ifacod[50] = 116 
ifacod[51] = 116 
ifacod[52] = 116 
ifacod[53] = 116 
ifacod[54] = 116 
#
#-- O3 HITRAN :
ifacod[55] = 20000 
ifacod[56] = 20000 
ifacod[57] = 20000 
#
#-- N2O HITRAN :
ifacod[58] = 20001 
ifacod[59] = 20001 
ifacod[60] = 20001 
ifacod[61] = 20001 
ifacod[62] = 20001 
#
#-- CO HITRAN :
ifacod[63] = 107 
ifacod[64] = 107 
ifacod[65] = 107 
ifacod[66] = 107 
ifacod[67] = 107 
#
#-- CH4 HITRAN :
ifacod[68] = 139 
ifacod[69] = 139 
ifacod[70] = 139 
#
#-- O2 HITRAN :
ifacod[71] = 117 
ifacod[72] = 117 
ifacod[73] = 117 
#-- NO HITRAN :
ifacod[74] = 115 
ifacod[75] = 115 
ifacod[76] = 115 
#
#-- SO2 HITRAN :
ifacod[77] = 20002 
ifacod[78] = 20002 
#
#-- NO2 Hitran :
ifacod[79] = 20003 
#
#-- NH3 HITRAN :
ifacod[80] = 170 
ifacod[81] = 170 
#
#-- HNO3 HITRAN :
ifacod[82] = 20004 
#
#-- OH HITRAN :
ifacod[83] = 102 
ifacod[84] = 102 
ifacod[85] = 102 
#
#-- HF HITRAN :
ifacod[86] = 124 
#
#-- HCl HITRAN :
ifacod[87] = 123 
ifacod[88] = 123 
#
#-- HBr HITRAN :
ifacod[89] = 20005 
ifacod[90] = 20005 
#
#-- HI HITRAN :
ifacod[91] = 20006 
#
#-- ClO HITRAN :
ifacod[92] = 20007 
ifacod[93] = 20007 
#
#-- OCS HITRAN :
ifacod[94] = 155 
ifacod[95] = 155 
ifacod[96] = 155 
ifacod[97] = 155 
#
#-- H2CO HITRAN :
ifacod[98] = 20008 
ifacod[99] = 20008 
ifacod[100] = 20008 
#
#-- HOCl HITRAN :
ifacod[101] = 20009 
ifacod[102] = 20009 
#-- N2 HITRAN :
ifacod[103] = 114 
#
#-- HCN HITRAN :
ifacod[104] = 137 
ifacod[105] = 137 
ifacod[106] = 137 
#
#-- CH3Cl HITRAN :
ifacod[107] = 20010 
ifacod[108] = 20010 
#
#-- H2O2 HITRAN :
ifacod[109] = 20011 
#
#-- C2H2 HITRAN :
ifacod[110] = 138 
ifacod[111] = 138 
#
#-- C2H6 HITRAN :
ifacod[112] = 20012 
#
#-- PH3 HITRAN :
ifacod[113] = 20013 
#
#-- COF2 HITRAN :
ifacod[114] = 20014 
#
#-- SF6 HITRAN :
ifacod[115] = 20015 
#
#-- H2S HITRAN :
ifacod[116] = 154 
ifacod[117] = 154 
ifacod[118] = 154 
#
#-- CO NASA-Goorvitch :
ifacod[119] = 107 
ifacod[120] = 107 
ifacod[121] = 107 
ifacod[122] = 107 
ifacod[123] = 107 
ifacod[124] = 107 
ifacod[125] = 107 
#
#-- CN Joergensen
ifacod[126] = 106 
ifacod[127] = 106 
#
# SiO
ifacod[128] = 121
ifacod[129] = 121    
ifacod[130] = 121    
#
# H3+
ifacod[131] = 21301  
#
# YO
ifacod[132] = 133    
#
# ZrO
ifacod[133] = 118    
ifacod[134] = 118    
ifacod[135] = 118    
ifacod[136] = 118    
ifacod[137] = 118    
ifacod[138] = 118    
ifacod[139] = 118    
#
# H2O [AMES]
ifacod[140] = 112    
#
# CH4 [GEISA]
ifacod[141] = 139    
ifacod[142] = 139    
ifacod[143] = 139    
ifacod[144] = 139    
ifacod[145] = 139    
ifacod[146] = 139    
#
#-- TiO [Schwenke]:
ifacod[147] = 111 
ifacod[148] = 111 
ifacod[149] = 111 
ifacod[150] = 111 
ifacod[151] = 111 

#
#-- FeH [Davis/Allard]:
ifacod[152] = 194 
#
#-- H16OD-AMES
ifacod[153] = 112 
#
#-- CH4 STDS [Champion & Hillico]
ifacod[154:161] = [139]*8

#
#-- VO AMES
ifacod[162] = 119 
#
#-- CrH AMES
ifacod[163] = 40600 

#
#-- H2O [SCAN]
ifacod[164] = 112 

#-- CH4 STDS199 [Champion & Hillico / Boudon]
#    
ifacod[165:186] = [139]*22


#
#-- MgH Phillip & Phillip
ifacod[187:189] = [108]*3

#
#-- CaH Phillip & Phillip
ifacod[190:191] = [109]*2 

#
#-- H2O [BT2]
ifacod[192] = 112 

#
#-- TiO [Plez2004]:
ifacod[193:197] = [111]*5

#
#-- VO [Plez2004]
ifacod[198:198] = [119]

#
#-- H2O [HITRAN2004]
ifacod[199:204] = [112]*6

#
#-- CO2 [HITRAN2004]
ifacod[205:212] = [116]*8

#
#-- O3 [HITRAN2004]
ifacod[213:217] = [59500]*5

#
#-- N2O [HITRAN2004]
ifacod[218:222] = [58800]*5

#
#-- CO [HITRAN2004]
ifacod[223:228] = [107]*6

#
#-- CH4 [HITRAN2004]
ifacod[229:231] = [139]*3

#
#-- O2 [HITRAN2004]
ifacod[232:234] = [117]*3

#
#-- NO [HITRAN2004]
ifacod[235:237] = [115]*3

#
#-- SO2 [HITRAN2004]
ifacod[238:239] = [29300]*2

#
#-- NO2 [Hitran2004]
ifacod[240:240] = [30000]

#
#-- NH3 [HITRAN2004]
ifacod[241:242] = [170]*2

#
#-- HNO3 [HITRAN2004]
ifacod[243:243] = [75200]

#
#-- OH [HITRAN2004]
ifacod[244:246] = [102]*3

#
#-- HF [HITRAN2004]
ifacod[247:247] = [124]

#
#-- HCl [HITRAN2004]
ifacod[248:249] = [123]*2

#
#-- HBr [HITRAN2004]
ifacod[250:251] = [36700]*2

#
#-- HI [HITRAN2004]
ifacod[252:252] = [44300]

#
#-- ClO [HITRAN2004]
ifacod[253:254] = [39700]*2

#
#-- OCS [HITRAN2004]
ifacod[255:259] = [155]*5

#
#-- CO2 [HITEMP]
ifacod[319:326] = [116]*8

#
del ifacod[327:]