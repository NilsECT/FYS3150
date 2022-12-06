from pyarma import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

A = cube()
A.load("test.bin", arma_binary)
print(size(A))
# A.reshape(198,320,198)

B = mat()
B.load("pot.bin", arma_binary)

# # print(A[:,:,0])
# P = np.zeros((198,198,320))
# for t in range(320):
#     for i in range(198):
#         for ii in range(198):
#             P[i,ii,t] = A[i,ii,t]

P = np.array(A)
V = np.array(B)
print(P.shape)
print(V.shape)

# for i in range(len(P[:,0,0])):
#     plt.figure()
#     plt.contourf(P[i])
#     plt.show()


#
# Let's generate a dummy time series for a function z(x,y,t)
#

# Set up a 2D xy grid

h = 0.005
x_points = np.arange(0, 1+h, h)
y_points = np.arange(0, 1+h, h)
x, y = np.meshgrid(x_points, y_points, sparse=True)

# plt.contourf(P[0])
# plt.colorbar()
# plt.legend()
# plt.show()

# Array of time points
dt = 0.005
t_points = np.arange(0, 1+dt, dt)



z_data_list = []
for t in range(320):
    z_data = P[t,:,:] + V/1e13
    z_data_list.append(z_data)

#
# Now the list z_data_list contains a series of "frames" of z(x,y,t), 
# where each frame can be plotted as a 2D image using imshow. Let's
# animate it!
#

# Some settings
fontsize = 12
t_min = t_points[0]
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]

# Create figure
fig = plt.figure()

ax = plt.gca()

# Create a colour scale normalization according to the max z value in the first frame
norm = mpl.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

# Plot the first frame
img = ax.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

# Axis labels
plt.xlabel("x", fontsize=fontsize)
plt.ylabel("y", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label("z(x,y,t)", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    norm = mpl.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))
    img.set_norm(norm)

    # Update z data
    img.set_data(z_data_list[i])

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=True, blit=0)

# Run the animation!
# plt.show()

# # Save the animation
anim.save('./animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)
