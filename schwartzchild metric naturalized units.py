import numpy as np
import matplotlib.pyplot as plt

M = 4.434*(10**(-3))#mass earth
M = 1.476*(10**(3))#mass sun
G = 1

pos = np.zeros(4)
pos[2] = 3.14159/2 #going about the equator
pos[1] = 400*1000 + 6.37*(10**6) #distance from center of earth to iss
pos[1] = 58*1000*(10**6)#distance from center of sun to mercury

vel = np.zeros(4)
vel[0] = 1
vel[3] = 0.0010895/(3*(10**8))#velocity of iss
vel[3] = 8.264*(10**(-7))/(3*(10**8))

acc = np.zeros(4)
index = ["t", "r", "h", "p"]

rs = 2*G*M

def explicit_der(current_pos, a = 0, b = 0, c = 0):
	return 0

def metric(current_pos, a, b):

	t = current_pos[0]
	r = current_pos[1]
	h = current_pos[2]
	p = current_pos[3]
	
	if a == 0 and b == 0:
		return -1 + rs/r
	elif a == 1 and b == 1:
		return (1 - rs/r)**(-1)
	elif a == 2 and b == 2:
		return r**2
	elif a == 3 and b == 3:
		return (r*np.sin(h))**2
	else:
		return 0

#change this to be the small modification of minkowski metric and see what happens, go be a gamer
def del_metric(current_pos, a = 0, b = 0, c = 0, d = 0.00001):
	#a and b are indices of the metric, c is the coordinate the derivative is taken wrt
	#d is the incriment to do derivative/integrals

	del_vec = np.zeros(4)
	del_vec[c] = d
	return (metric(current_pos, a, b) - metric(current_pos+del_vec, a, b))/d

def inv_metric(current_pos):
	metric_at_point = np.zeros((4, 4))
	for a in range(4):
		for b in range(4):
			metric_at_point[a][b] = metric(current_pos, a, b)
	#print(metric_at_point)
	return np.linalg.inv(metric_at_point)

def calc_chris(current_pos, g_inv, a, b, c):

	chris_cab = 0

	for d in range(4):
		g_cd = g_inv[c][d]
		if d > 1:
			delta = 0.001
		else:
			delta = 10000
		chris_cab += 0.5*g_cd*(del_metric(current_pos, d, b, a, delta) + del_metric(current_pos, d, a, b, delta) - del_metric(current_pos, a, b, d, delta))

	return chris_cab

def find_acc(current_pos, current_vel):
	del_acc = np.zeros(4)
	for c in range(4):
		g_inv = inv_metric(current_pos)
		for a in range(4):
			for b in range(4):
				del_acc[c] += 1*calc_chris(current_pos, g_inv, a, b, c)*current_vel[a]*(current_vel[b])
				#some manner of extreme buffoonery has come over me and I am at a loss for words
				#the equation above should have a negative sign in it
				#and yet, when it does, the model does not generate physically realistic results
				#which means I have missed another negative sign in another area
				#I will fix this later
				
	return del_acc


xs = []
ts = []
print("started")
for t in range(24*365*100):
	dt = 60*60
	acc = find_acc(pos, vel)
	k1v = acc*dt*(3*10**8)
	k1x = vel*dt*(3*10**8)
	k2v = find_acc(pos + k1x/2, vel + k1v/2)*dt*(3*10**8)
	k2x = (vel + k1v/2)*dt*(3*10**8)
	pos += k2x
	vel += k2v

	#print(vel)
	if (t%1000 == 0):
		print(t)
	#xs.append(acc[1])
	#ts.append(t)
	xs.append(pos[1]*np.sin(pos[3]))
	ts.append(pos[1]*np.cos(pos[3]))
plt.plot(ts, xs)
plt.show()




