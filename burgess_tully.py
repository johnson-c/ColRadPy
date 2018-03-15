#dipole transitions
#type 1
dp_inds = np.where(ne['inf_engy'] <0)[0]

c = 1.5
x = 1 - np.log(c) / np.log( ne['temp_grid'][-1]*0.69488623/ (ne['energy'][ne['col_transitions'][:,0]-1] - ne['energy'][ne['col_transitions'][:,1] -1]) + c)
y = ne['col_excit'][:,-1] / np.log( ne['temp_grid'][-1]*0.69488623/ (ne['energy'][ne['col_transitions'][:,0]-1] - ne['energy'][ne['col_transitions'][:,1] -1])  + np.exp(1))

y_1 = np.abs( ne['inf_engy'][dp_inds])

m = (y_1 - y[dp_inds]) / (1 -x[dp_inds])
slope = m
b =  y[dp_inds] - slope*x[dp_inds]

x_arr_dp = 1- np.log(c) / np.log( ne['temp_grid'].reshape(len(ne['temp_grid']),1)*0.69488623/ (ne['energy'][ne['col_transitions'][:,0]-1] - ne['energy'][ne['col_transitions'][:,1] -1]) + c)

#ne['col_excit'][:,-1] changed 3/11/17
y_arr_dp = ne['col_excit'][:,:] / np.log( ne['temp_grid'].reshape(len(ne['temp_grid']),1)*0.69488623/ (ne['energy'][ne['col_transitions'][:,0]-1] - ne['energy'][ne['col_transitions'][:,1] -1])  + np.exp(1))



for i in range(0,len(dp_inds)):
    plt.figure()
    plt.title('Transition ' +str(dp_inds[i]+1) + ',  ' + str(ne['col_transitions'][dp_inds[i],0]) + '  ->  ' +str(ne['col_transitions'][dp_inds[i],1]))
    plt.plot(x_arr_dp[:,dp_inds[i]],y_arr_dp[:,dp_inds[i]],color='b',label='Data')
    plt.scatter( np.array([x[dp_inds][i],1]), np.array([y[dp_inds][i],y_1[i]]),label='Fit points')
    plt.plot( np.linspace(x[dp_inds][i],1,100), slope[i] * np.linspace(x[dp_inds][i],1,100) +b[i],color='r',label='Infinite Energy fit')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend(loc=2)
    plt.savefig('/home/curtis/lewis/ne_1/dp_transition_'+ str(dp_inds[i]+1) + '.png')
    print(i)
