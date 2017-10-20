PRO inspect,field,object,setup,mode,ps=ps,norm=norm,dir=dir,yr=yrp,labels=labels

if not keyword_set(dir) then spec='SPECTRA/' else spec='SPECTRA_'+dir+'/'
if not keyword_set(dir) then dir='OUTPUT/' else dir='OUTPUT_'+dir+'/'

dum10 = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(dum10)/2., SIN(dum10)/2., /FILL
colors

file=dir+field+'_'+object+'_'+setup+'_'+mode+'_SME.out'

if not file_test(file) then begin 
   print,dir+field+'_'+object+'_'+setup+'_'+mode+'_SME.out not found'
   return
endif

restore,file
if keyword_set(norm) and mode eq 'Sp' then begin 
   file=findfile(spec+field+'_'+object+'_'+setup+'_Sp.dat')
   if file[0] ne '' then begin 
      readcol,file[0],wave_norm,x,x,smod_norm,format='d,d,d,d' ,/silent 
   endif else begin 
      file=findfile(spec+field+'_'+object+'_'+setup+'.dat')
      if file[0] ne '' then readcol,file[0],wave_norm,x,x,smod_norm,format='d,d,d,d' ,/silent 
   endelse
endif
if keyword_set(norm) and mode ne 'Sp' then begin 
   file=findfile(spec+field+'_'+object+'_'+setup+'_'+mode+'.dat')
   if file[0] ne '' then readcol,file[0],wave_norm,x,x,smod_norm,format='d,d,d,d',/silent 
endif

if keyword_set(ps) then begin
   psfull,'l',name=field+'_'+object+'_'+setup+'_'+mode+'.ps'
   !p.multi=[0,1,2,0,0]
   maskcol=5
   !p.thick=3
endif else begin
   window,xsize=1200,ysize=400 
   maskcol=3
   !p.thick=2
endelse
   
marcs_abund,solar,fehmod=0.

elstr=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cs','Es']

for i=0,sme.nseg-1 do begin

   j=where(sme.wave ge sme.wran[0,i] and sme.wave le sme.wran[1,i] and sme.sob gt 0.,jc)

   if not keyword_set(yrp) then yr=minmax([sme.sob[j],sme.smod[j]]) else yr=yrp

   if keyword_set(wave_norm) and not keyword_set(yrp) then begin
      k=where(wave_norm ge sme.wran[0,i] and wave_norm le sme.wran[1,i])
      yr=minmax([sme.sob[j],sme.smod[j],smod_norm[k]])
   endif

   dy=[yr[1]-yr[0]]*0.05

   if sme.maxiter eq 20 then conv='' else conv=' NOT CONVERGED'

   delt=sme.wave[1:*]-sme.wave[0:n_elements(sme.wave)-2]
   delt=delt[where(sme.wave ge sme.wran[0,i] and sme.wave le sme.wran[1,i])]
   frees = where(sme.ab_free,frees_i)
   ;print,'Free ab: ',frees
   e1 = '  '
   e2 = '  '
   e3 = '  '
   abf1='-    '
   abf2='-    '
   abf3='-    '
   if max(frees) ne -1 then begin
      if frees_i ge 1 then abf1=string(sme.abund[frees[0]]-solar[frees[0]],format='(f5.2)')
      if frees_i ge 1 then e1 = elstr[frees[0]]
      if frees_i ge 2 then abf2=string(sme.abund[frees[1]]-solar[frees[1]],format='(f5.2)')
      if frees_i ge 2 then e2 =elstr[frees[1]]
      if frees_i ge 3 then abf3=string(sme.abund[frees[2]]-solar[frees[2]],format='(f5.2)')
      if frees_i ge 3 then e3 =elstr[frees[2]]
   endif
   line=where(sme.mob eq 1 and sme.wave ge sme.wran[0,i] and sme.wave le sme.wran[1,i],linec)
   cont=where(sme.mob eq 2 and sme.wave ge sme.wran[0,i] and sme.wave le sme.wran[1,i],contc)
   segm=where(sme.wran[0,i] and sme.wave le sme.wran[1,i],segmc)
   nan=where(sme.wave[1:*]-sme.wave[0:*] ge min(delt)*1.5 and sme.wave ge sme.wran[0,i] and sme.wave le sme.wran[1,i],nanc)

   sn=mean(sme.sob[line]/sme.uob[line])
   ;print,sn
   if mode ne 'Sp' then sn=min(sme.smod[line])*sn
   dy=[yr[1]-yr[0]]*0.049
   plot,[-99],xr=sme.wran[*,i],xs=1,yr=yr+[-dy,dy],ys=1,title=object+string(sme.teff,sme.grav, sme.feh, sme.vmic,e1,abf1,e2,abf2,e3,abf3,sn,format='(": T=",i4, " K, g=", F5.2, ", m=", F+5.2,", x=", F+5.2,", ", a2, "=", a5, ", ", a2, "=", a5, ", ",a2,"=", a5, ", d=", F3.0)')+conv


   if linec ne 0 then begin
      for k=0,linec-2 do begin 
         x=sme.wave[line([k,k+1,k+1,k,k])]
         y=[yr[0]-dy,yr[0]-dy,yr[1]+dy,yr[1]+dy,yr[0]-dy]
         if line[k+1]-line[k] eq 1 then polyfill,x,y,col=maskcol
      endfor
   endif
   if contc ne 0 then begin
      for k=0,contc-2 do begin 
         x=sme.wave[cont([k,k+1,k+1,k,k])]
         y=[yr[0]-dy,yr[0]-dy,yr[1]+dy,yr[1]+dy,yr[0]-dy]
         if cont[k+1]-cont[k] eq 1 then polyfill,x,y,col=4
      endfor
   endif

   sortwave=sort(sme.wave)
   oplot,sme.wave[sortwave],sme.sob[sortwave],psym=8
   errplot,sme.wave[sortwave],sme.sob[sortwave]-sme.uob[sortwave],sme.sob[sortwave]+sme.uob[sortwave]
   oplot,sme.wave[sortwave],sme.smod[sortwave],col=2
   if keyword_set(wave_norm) then oplot,wave_norm,smod_norm,col=4

   if keyword_set(labels) then begin 
      if labels eq 1 then labels=9
      j=where(sme.atomic[2,*] gt sme.wran[0,i] and sme.atomic[2,*] lt sme.wran[1,i],jc)
      j=j[reverse(sort(sme.depth[j]))]
      j=j[0:min([labels-1,jc-1])]
      j=j[sort(sme.atomic[2,j])]
      for jj=0,n_elements(j)-1 do begin 
         x=sme.atomic[2,j[jj]]
         ;print,x,sme.species[j[jj]]
         y=1-sme.depth[j[jj]]
         oplot,[x,x],[y,y-0.05]
         xyouts,x,y-0.15,sme.species[j[jj]],orientation=90,charsize=1.0
      endfor
   endif

   if not keyword_set(ps) then read,ok

endfor

!p.multi=0
if keyword_set(ps) then psoff

END
