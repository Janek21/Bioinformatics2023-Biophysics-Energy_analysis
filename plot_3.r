setwd('./Biophysics_A1')
table <- read.csv('inter_en_res.csv')

#plot of electrostatic
combined_1<-(cbind(table$elec_res, table$elec_ala))
colnames(combined_1) <- c('Electrostatic of Result', 'Electrostatic with Alanine')
ts.plot(combined_1, col=c("blue", "orange"), ylab='Value', xlab='Residue', main = 'Electrostatic Interactions' )
legend("topright", legend = colnames(combined_1), col = c('blue', 'orange'), lty = 1,)

#plot of vdw
combined_2<-(cbind(table$vdw_res, table$vdw_ala))
colnames(combined_2) <- c('Van der Waals of Result', 'Van der Waals with Alanine')
ts.plot(combined_2, col=c("blue", "orange"), ylab='Value', xlab='Residue', main = 'Van der Waals Interactions')
legend("bottomright", legend = colnames(combined_2), col = c('blue', 'orange'), lty = 1,)

#plot of solvation of A-B interaction
combined_3<-(cbind(table$solv_AB_res, table$solv_AB_ala))
colnames(combined_3) <- c('Solvation of A-B interaction of Result', 'Solvation of A-B interaction with Alanine')
ts.plot(combined_3, col=c("blue", "orange"), ylab='Value', xlab='Residue', main = 'Solvation A-B interaction')
legend("topleft", legend = colnames(combined_3), col = c('blue', 'orange'), lty = 1,)

#plot of solvation of A
combined_4<-(cbind(table$solv_A_res, table$solv_A_ala))
colnames(combined_4) <- c('Solvation of A of result', 'Solvation of A with Alanine')
ts.plot(combined_4, col=c("blue", "orange"), ylab='Value', xlab='Residue', main = 'Solvation of A')
legend("topleft", legend = colnames(combined_4), col = c('blue', 'orange'), lty = 1,)
