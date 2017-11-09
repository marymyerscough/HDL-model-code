name = 'athst1wHDL'

#r3 = run(e=name,c=name,DS='-',UZSTOP={'int(N)':[0,1],'Sigma_h':[0,800]})
#filename = 'alpha_m_1e-1'
#r3.writeRawFilename(filename)

#mv alpha_m_1e-1 Lo_l_Hi_m

r1 = run(e=name,c=name,DS='-',STOP=['LP1'])
r11 = run(r1('LP1'),ICP=['Sigma_h','alpha_m'],UZSTOP={'alpha_m':[1.0e-1,1.0e2]},ISW=2)
r12 = run(r11,DS='-',DSMAX=1.5)
r2 = run(r12,ICP=['Sigma_h','alpha_l'],UZSTOP={'alpha_l':[1.0e-1,1.0e1]},ISW=2,DSMAX=3.0)
#r3 = run(r2,DS='-',DSMAX=0.8)#
filename = 'test123'
r2.writeRawFilename(filename)

#mv p2_alpha_l=1.4e-1_sigma_h_alpha_m Hi_l_Hi_m