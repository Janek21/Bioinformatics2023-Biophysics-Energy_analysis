table <- read.csv('inter_en_res.csv')


combined_1<-(cbind(table$elec_res, table$elec_ala))
ts.plot(combined, col=c("blue", "orange"), ylab='Value', xlab='Residue', main = 'Electrostatic Interactions')

combined_2<-(cbind(table$vdw_res, table$vdw_ala))
ts.plot(combined_2, col=c("blue", "orange"), ylab='Value', xlab='Residue', main = 'Van der Waals Interactions')

combined_3<-(cbind(table$solv_AB_res, table$solv_AB_ala))
ts.plot(combined_3, col=c("blue", "orange"), ylab='Value', xlab='Residue', main = 'Solvation A-B interaction')

combined_4<-(cbind(table$solv_A_res, table$solv_A_ala))
ts.plot(combined_4, col=c("blue", "orange"), ylab='Value', xlab='Residue', main = 'Solvation of A')
