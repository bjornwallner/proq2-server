set terminal png font times
set output '/home/bjowa51/proq2-server/server/output/2017-2/job.sQsvy6q6/svmdata-Sscore3-ac4-atom1-entropy3-grsa_sc1-gss_sc1-profile0-prsa0-pw1-pwin23-rc6-res1-restypes6-rsa13-rsa_sc21-ss1-ss_sc21-stride5-surf1001-surf251-surf501-surf751-termini5-topology0-z0-z_sc0-zpred0/input.pdb.svm.pred.png'
set xlabel 'sequence position'
set ylabel 'predicted quality, 1-best 0-worst'
plot '/home/bjowa51/proq2-server/server/output/2017-2/job.sQsvy6q6/svmdata-Sscore3-ac4-atom1-entropy3-grsa_sc1-gss_sc1-profile0-prsa0-pw1-pwin23-rc6-res1-restypes6-rsa13-rsa_sc21-ss1-ss_sc21-stride5-surf1001-surf251-surf501-surf751-termini5-topology0-z0-z_sc0-zpred0/input.pdb.svm.pred.all' u 1:2 w l lw 1 title 'raw','' u 1:2 w l lw 2 smooth bezier title 'smoothed'
set ylabel 'predicted distance deviation (angstroms)'
set output '/home/bjowa51/proq2-server/server/output/2017-2/job.sQsvy6q6/svmdata-Sscore3-ac4-atom1-entropy3-grsa_sc1-gss_sc1-profile0-prsa0-pw1-pwin23-rc6-res1-restypes6-rsa13-rsa_sc21-ss1-ss_sc21-stride5-surf1001-surf251-surf501-surf751-termini5-topology0-z0-z_sc0-zpred0/input.pdb.svm.pred.D.png'
plot '/home/bjowa51/proq2-server/server/output/2017-2/job.sQsvy6q6/svmdata-Sscore3-ac4-atom1-entropy3-grsa_sc1-gss_sc1-profile0-prsa0-pw1-pwin23-rc6-res1-restypes6-rsa13-rsa_sc21-ss1-ss_sc21-stride5-surf1001-surf251-surf501-surf751-termini5-topology0-z0-z_sc0-zpred0/input.pdb.svm.pred.all' u 1:3 w l lw 1 title 'raw','' u 1:3 w l lw 2 smooth bezier title 'smoothed'
