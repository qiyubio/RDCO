norfiles=dir('simulation_canfam3.bed_chr*.bed')
for file =norfiles'
     cal_APCA_new(file.name)
end
