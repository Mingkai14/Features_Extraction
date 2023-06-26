library("bio3d")

arg=commandArgs(T)
wt_pdb_path<-arg[1]
mut_pdb_path<-arg[2]
loc<-arg[3]
path<-arg[4]
outpath<-arg[5]

if (is.na(wt_pdb_path)|is.na(mut_pdb_path)|is.na(loc)|is.na(path)|is.na(outpath))
{
   print('Error!')
   quit()
}

#wt_pdb<-read.pdb(wt_pdb_path)
#wt_modes <- nma(wt_pdb)
#plot(wt_modes)
#print(wt_modes$fluctuations)
#print(wt_modes$fluctuations[1])
#print(wt_modes$fluctuations[42])
#print(length(wt_modes$fluctuations))
#wt_energy<-deformation.nma(wt_modes)
#print(wt_energy$ei[42,])
#print(wt_energy$ei[1,])
#print(wt_energy$ei[2,])
#print(wt_energy$sums)


#mut_pdb<-read.pdb(mut_pdb_path)
#mut_modes <- nma(mut_pdb)
#plot(mut_modes)
#print(mut_modes$fluctuations)
#print(mut_modes$fluctuations[1])
#print(mut_modes$fluctuations[42])
#print(length(mut_modes$fluctuations))
#mut_energy<-deformation.nma(mut_modes)
#print(mut_energy$ei[42,])
#print(mut_energy$ei[1,])
#print(mut_energy$ei[2,])
#print(mut_energy$sums)


#r <- rmsip(wt_modes,mut_modes,subset=10, row.name="a", col.name="b")


wt_pdb<-read.pdb(wt_pdb_path)
mut_pdb<-read.pdb(mut_pdb_path)
l<-list(wt_pdb,mut_pdb)
pdbs <- pdbaln(l,exefile= paste(path, "muscle", sep=""),outfile=paste(outpath,"aln.fa",sep=""))
modes <- nma(pdbs, rm.gaps=TRUE)
loc_<-as.numeric(loc)
f_wt_loc<-as.double(modes$fluctuations[,loc_][1])
f_mut_loc<-as.double(modes$fluctuations[,loc_][2])
rmsip<-as.double(modes$rmsip[,1][2])
path="r_output.txt"
cat(f_wt_loc,f_mut_loc,rmsip,file=paste(outpath,path,sep=""))










