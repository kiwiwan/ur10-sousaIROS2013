# Initializations
datafolder = 'data/'
tmpfolder = 'tmp/'

from sympy import init_printing
init_printing()


import numpy

proc = numpy.load(tmpfolder + 'procdata/proc.npz')
dq = proc['dq']
regr = numpy.load(tmpfolder + 'procdata/regr.npz')
omega = numpy.matrix(regr['omega'])
W = numpy.matrix(regr['W'])
del proc, regr


w = omega
v = dq.flatten()
s = len(v)


def selected_regression( obs_mat, output_vec, select_list=None ):
    
    if select_list is None:
        A = obs_mat[:,:]
        b = output_vec[:,:]
    else:
        A = obs_mat[select_list,:]
        b = output_vec[select_list,:]
    
    num = len(b)
    q1,r1 = numpy.linalg.qr(A)
    b1 = q1.T*b
    cond = numpy.linalg.cond(r1)
    x = numpy.linalg.inv(r1.T * r1) * r1.T * b1
    se = b - A*x
    err = numpy.linalg.norm(se)
    relerr = err / numpy.linalg.norm(b)
    mse = err**2/(num-A.shape[1])
    
    return x,se,{'num':num,'cond':cond,'mse':mse, 'relerr':relerr}


list_v_threshold = numpy.arange(0.0, 0.02, 0.002).tolist() + numpy.arange(0.02, 0.08+1e-10, 0.02).tolist()

selections = []
list_props = []
v_db_unc_regr = []
v_num = []
v_cond = []
v_relerr = []
v_mse = []
for idx, v_threshold in enumerate(list_v_threshold):
    r_vel = []
    if True:
        for i,vi in enumerate(v):
            if numpy.abs(vi) >= v_threshold:
                r_vel.append( i )
    else:                
        for i in range(len(v)/rbt.dof):
            elimin = False
            for l in range(rbt.dof):
                vi = v[i*rbt.dof+l]
                if numpy.abs(vi) < v_threshold:
                    elimin = True
            if not elimin:
                r_vel += [i*rbt.dof+l for l in range(rbt.dof)]
    selections.append(r_vel)
    db_unc_regr_vel,se,props = selected_regression( W, w, r_vel )
    list_props.append(props)
    v_db_unc_regr.append(db_unc_regr_vel)
    v_num.append(props['num'])
    v_cond.append(props['cond'])
    v_relerr.append(props['relerr'])
    v_mse.append(props['mse'])
    
    print "%d %.3f %s"%(idx, v_threshold, props)



from matplotlib import pyplot as plt

plt.figure(figsize=(5,3))
plt.plot(list_v_threshold,v_cond/v_cond[0]*100, label='Condition number')
plt.plot(list_v_threshold,v_mse/v_mse[0]*100, label='MSE')
#plt.plot(list_v_threshold,v_relerr/v_relerr[0]*100, label='rel err')
plt.legend(loc='right')
plt.axis([0.0,0.08,70,120])
plt.ylabel('Ratio to values of complete data set (%)')
plt.xlabel('Velocity threshold ($\mathrm{rad}/s$)')

plt.show()

threshold_idx = 5

(v_cond[threshold_idx] / v_cond[0]  - 1) * 100.0

(v_mse[threshold_idx]/v_mse[0] - 1) *100

v_cond[threshold_idx]

select = selections[threshold_idx]
threshold = list_v_threshold[threshold_idx]

W_n0v = W[select,:]
omega_n0v = omega[select,:]
    
Q1_n0v,R1_n0v = numpy.linalg.qr(W_n0v)
rho1_n0v = Q1_n0v.T*omega_n0v

numpy.savez_compressed(tmpfolder + 'procdata/regr_n0v', W=W_n0v, omega=omega_n0v, Q1=Q1_n0v, R1=R1_n0v, rho1=rho1_n0v)

del W, omega, W_n0v, omega_n0v, Q1_n0v, R1_n0v, rho1_n0v


