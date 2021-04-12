
importCommon()

random_list = []
for i in range(35):
    if i<7:
        random_list.append( numpy.random.normal(0.6, 0.1) )
    elif i<14:
        random_list.append( numpy.random.normal(0.2, 0.1) )
    elif i<21:
        random_list.append( numpy.random.normal(0.9, 0.3) )
    elif i<28:
        random_list.append( numpy.random.normal(0.2, 0.1) )
    else:
        random_list.append( numpy.random.normal(0.6, 0.1) )

random_list = numpy.array(random_list)

colors = numpy.array( ['#000000']*35, dtype=str )
colors[ random_list>0.4 ] = "#ff9800" 
colors[ random_list>0.7 ] = "#ff0000" 

plt.bar(range(len(random_list)), random_list, color=colors)
plt.savefig("/tmp/fig.pdf")
plt.show()


NAI = []
for i in range(35):
    if i<7:
        NAI.append( numpy.random.normal(0.6, 0.1) )
    elif i<14:
        NAI.append( numpy.random.normal(0.2, 0.1) )
    elif i<21:
        NAI.append( numpy.random.normal(0.9, 0.3) )
    elif i<28:
        NAI.append( numpy.random.normal(0.2, 0.1) )
    else:
        NAI.append( numpy.random.normal(0.6, 0.1) )

DMSO = []
for i in range(35):
    DMSO.append( numpy.random.normal(0.2, 0.05) )

plt.step(range(len(NAI)), NAI, color='red')
plt.step(range(len(DMSO)), DMSO, color='gray')

plt.savefig("/tmp/fig.pdf")
plt.show()
