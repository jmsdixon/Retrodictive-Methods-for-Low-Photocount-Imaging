# Retrodictive-Methods-for-Low-Photocount-Imaging
This package contains three low photocount image reconstruction methods.  They comprise a mix of smoothing, quantum retrodiction and machine learning. One is a replication of work by Dr Grame Johnstone and the files '', '' and '' are his work.  The rest is entirely my own work.

This package contains 3 low intensity light imaging algorithms.  2 methods I developed myself and the 3rd was produced by Dr. Graeme Johnstone of the Institute of Photonics at the University of Strathclyde and is a replication of the work of Speiret’s et. al. (2017) “From Retrodiction to Bayesian Quantum Imaging”. 
The files here can be used to simulate a low intensity light imaging of a PNG of the user’s choice, set simulate_t = 1 and a PNG in Interface.m and run.
The files Bayret_mat_Mj.m, Bayret.m were developed by Dr Johnstone and shared with me – they are uploaded with his permission.  The file Bayret_mat_Mj_fast.m is originally Dr Johnstone’s but with a minor adaption by me.  The rest is entirely my own work.
A user will require their own image to test, I have provided my training data for individual frame, faces_ind.m, meaning a maximum of 2 counts per pixel per frame, requiring fpf=2; and effective frames of size 20 frames, faces_fsz20.m.  From these, my results for my 5th year Physics project were produced using the lab data provided by Dr Johnstone and test images.  For access to these images and the lab data please request it directly from myself.
Any questions or problems running the code, fire me an email I’m happy to assist.
