import numpy as np
from matplotlib import pyplot as plt

# converts degrees to radians
def deg_to_rad(deg):
    return (np.pi/180)*deg


# converts radians to degrees
def rad_to_deg(rad):
    return (180/np.pi)*rad


# converts a string input of hms ('00:00:00') into a float output of degrees (0)
def hms_to_deg(hms):

# splits the input by colons into its 3 sections
    hms_split = hms.split(':')

# converts from string to float and hms to deg, taking into account negative values
    if '-' in hms_split[0]:
        deg = (float(hms_split[0]) - (float(hms_split[1])/60) - (float(hms_split[2])/3600))*15
    else:
        deg = (float(hms_split[0]) + (float(hms_split[1])/60) + (float(hms_split[2])/3600))*15

# returns degree value
    return deg


# converts a string input of dms ('00:00:00') into a float output of degrees (0)
def dms_to_deg(dms):

# splits the input by colons into its 3 sections
    dms_split = dms.split(':')

# converts from string to float and dms to deg, taking into account negative values
    if '-' in dms_split[0]:
        deg = (float(dms_split[0]) - (float(dms_split[1])/60) - (float(dms_split[2])/3600))
    else:
        deg = (float(dms_split[0]) + (float(dms_split[1])/60) + (float(dms_split[2])/3600))

# returns degree value
    return deg


'''
tests:
1. whether or not the difference between the expected and actual values is less than the absolute tolerance
2. whether or not the expected value is equal to 0
3. whether or not the fractional difference is less than the relative tolerance
returns a boolean value
'''
def approx_equals(expected: float, actual: float, rel_tol: float = 1e-10, abs_tol: float = 1e-10):
    if abs(expected - actual) < abs_tol:
        return True
    elif expected == 0:
        return False
    elif abs((actual - expected)/expected) < rel_tol:
        return True
    else: 
        return False


'''
test:
1. whether or not the shapes of the 2 arrays are equal
2. whether or not each item in each array is approximately equal to one another based on "approx_equals" function
returns a boolean value
'''
def array_approx_equals(arr_expected: np.array, arr_actual: np.array, rel_tol: float = 1e-10, abs_tol = 1e-10):

    # 1.
    if np.shape(arr_expected) != np.shape(arr_actual):
        return False
    arr_expected.flatten()
    arr_actual.flatten()

    # 2. 
    for i in arr_expected:
        if not approx_equals(arr_expected[i], arr_actual[i], rel_tol, abs_tol):
            return False
        else:
            return True


# converts spherical coordinates to rectangular coordinates
def spherical_to_rectangular(spherical_coords: np.array):

    # converts each coordinate from degrees to radians and assigns it to a variable
    rho = deg_to_rad(spherical_coords[0])
    phi = deg_to_rad(spherical_coords[1])
    theta = deg_to_rad(spherical_coords[2])

    # performs conversion from spherical coordinates to rectangular coordinates
    x = rho*np.sin(theta)*np.cos(phi)
    y = rho*np.sin(theta)*np.sin(phi)
    z = rho*np.cos(theta)

    # creates new array to store rectangular coordinates
    rect_coords = np.array([x,y,z])

    return rect_coords


# applies "spherical_to_rectangular" function to every set of spherical coordinates in a given Nx3 array
def spherical_to_rectangularNx3(spherical_coords_array):

    # prepares an empty list to append arrays of rectangular coordinates to
    rect_coords_list = []

    # iterates through spherical coordinates in spherical_coords_array, converts them to rectangular coordinates, 
    # and appends them to rect_coords_list
    for spherical_coords in spherical_coords_array:
        rect_coords_list.append(spherical_to_rectangular(spherical_coords))

    # returns rect_coords_list as an array
    return np.array(rect_coords_list)


# converts rectangular coordinates to spherical coordinates (inverse function of "spherical_to_rectangular")
def rectangular_to_spherical(rect_coords: np.array):

    # assigns each coordinate to a variable
    x = rect_coords[0]
    y = rect_coords[1]
    z = rect_coords[2]

    # performs conversion from rectangular coordinates to spherical coordinates
    rho = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    phi = rad_to_deg(y/z)
    theta = rad_to_deg(z/rho)

    # creates new array to store spherical coordinates
    spherical_coords = np.array([rho, phi, theta])

    return spherical_coords


# applies "rectangular_to_spherical" function to every set of rectangular coordinates in a given Nx3 array
def rectangular_to_sphericalNx3(rect_coords_array):

    # prepares an empty list to append arrays of spherical coordinates to
    spherical_coords_list = []

    # iterates through rectangular coordinates in rect_coords_array, converts them to spherical coordinates,
    # and appends them to spherical_coords_list

    for rect_coords in rect_coords_array:
        spherical_coords_list.append(rectangular_to_spherical(rect_coords))

    # returns spherical_coords_list as an array
    return np.array(spherical_coords_list)


# calculates LST using Julian date and longitude
def JD_to_LST(JD, long):

    # calculates JD number
    JD0 = int(JD-0.5) + 0.5

    # calculates day
    d = JD - JD0

    # calculates Julian centuries
    JD100 = ((JD0 - 2451545)/36525)

    # calculates GMST
    theta0 = (100.46061837 + 36000.77053608*JD100 - (3.87933*(10**-4)*(JD100**2)) - ((JD100*3)/(3.871*10**7)))
    thetaG = theta0 + 360.985647366*d

    # calculates LST
    LST = float(thetaG + long)%360

    return LST


# calculates azimuth and altitude using LST, RA, dec, and latitude
def azalt(LST_deg, RA_deg, dec_deg, lat_deg):

    # converts each value from degrees to radians
    LST = deg_to_rad(LST_deg)
    RA = deg_to_rad(RA_deg)
    dec = deg_to_rad(dec_deg)
    lat = deg_to_rad(lat_deg)
    
    # calculates hour angle
    H = LST - RA
    
    # calculates altitude
    h = np.arcsin(((np.sin(dec)*np.sin(lat))) + ((np.cos(dec))*(np.cos(lat))*(np.cos(H))))
    
    # calculates arcsin and arccos of azimuth and converts them from radians to degrees
    A_arcsin = rad_to_deg(np.arcsin((-np.cos(dec))*np.sin(H))/(np.cos(h)))
    A_arccos = rad_to_deg(np.arccos(np.sin(dec)-((np.sin(h))*np.sin(lat)))/(np.cos(h)*np.cos(lat)))
    
    # performs quadrant check to determine true value of azimuth
    if A_arccos > 90:
        A = 180 - A_arcsin
    elif A_arcsin < 0:
        A = A_arcsin
    else:
        A = 360 + A_arcsin
    
    # result is a tuple containing azimuth and altitude in degrees
    result = (A, rad_to_deg(h))

    return result


# calculates azimuth and altitude using JD, long, RA, dec, and lat (combines "JD_to_LST" and "azalt")
def JD_LST_azalt(JD, long, RA_deg, dec_deg, lat_deg):

    # calculates LST in degrees using "JD_to_LST" function
    LST_deg = (JD_to_LST(JD, long))

    # converts values from degrees to radians
    LST = deg_to_rad(LST_deg)
    RA = deg_to_rad(RA_deg)
    dec = deg_to_rad(dec_deg)
    lat = deg_to_rad(lat_deg)

    # calculates hour angle
    H = LST - RA

    # calculates altitude
    h = np.arcsin(((np.sin(dec)*np.sin(lat))) + ((np.cos(dec))*(np.cos(lat))*(np.cos(H))))

    # calculates arcsin and arccos of azimuth and converts them from radians to degrees
    A_arcsin = rad_to_deg(np.arcsin((-np.cos(dec))*np.sin(H))/(np.cos(h)))
    A_arccos = rad_to_deg(np.arccos(np.sin(dec)-((np.sin(h))*np.sin(lat)))/(np.cos(h)*np.cos(lat)))

    # performs quadrant check to determine true value of azimuth
    if A_arccos > 90:
        A = 180 - A_arcsin
    elif A_arcsin < 0:
        A = A_arcsin
    else:
        A = 360 + A_arcsin

    # result is a tuple containing azimuth and altitude in degrees
    result = (A, rad_to_deg(h))

    return result

    '''
    Alternate approach, using "azalt" function:
    LST_deg = (JD_to_LST(JD, long))
    result = azalt(LST_deg, RA_deg, dec_deg, lat_deg)
    return result
    '''


# calculates angular momentum using position and velocity vectors
def ang_momentum(pos, vel):

    # iterates through position and velocity vectors and converts values into AU and AU/Gaussian Day
    for i in range(len(pos)):
        pos[i] = pos[i]*(6.6845871222684*(10**-9))
    for i in range(len(vel)):
        vel[i] = vel[i]*(6.6845871222684*(10**-9))*(86400)*(58.13244087)
    
    # calculates angular momentum using cross product of position and velocity vectors
    ang_momentum = np.cross(pos,vel)

    # iterates through angular momentum vector and rounds values to 6 digits
    for i in range(len(ang_momentum)):
        ang_momentum[i] = round(ang_momentum[i], 6)
    
    return ang_momentum


# given a text file appropriately formatted for the "orbital_elements" function, extracts data
def oe_data(file):

# opens text file and creates a list of lines
    data = open(file)
    content = data.readlines()

# selects position vector (first item in the list) and splits it into its components
    r_content = content[0]
    r_split = r_content.split(' ')

# prepares an empty array to append position components to
    r = np.array([])

# iterates through position vector and converts each component into floats
    for i in range(len(r_split)):
        r_split[i] = float(r_split[i])
        r = np.append(r, r_split[i])

# selects velocity vector (second item in the list) and splits it into its componenets
    r_dot_content = content[1]
    r_dot_split = r_dot_content.split(' ')

# prepares an empty array to append velocity components to
    r_dot = np.array([])

# iterates through velocity vector and converts each component into floats
    for i in range(len(r_dot_split)):
        r_dot_split[i] = float(r_dot_split[i])
        r_dot = np.append(r_dot, r_dot_split[i])

# extracts and organizes JPL Horizons data
    JD = float(content[2])
    a_jpl = float(content[3])
    e_jpl = float(content[4])
    i_jpl = float(content[5])
    O_jpl = float(content[4])
    w_jpl = float(content[5])
    M_jpl = float(content[6])

# returns a list of data
    return [r, r_dot, JD, a_jpl, e_jpl, i_jpl, O_jpl, w_jpl, M_jpl]


'''
calculates orbital elements using position (r) and velocity (r_dot) vectors and compares them to JPL Horizons data
a: semimajor axis 
e: eccentricity
i: inclination
O: right ascenscion of the ascending node
w: argument of the perihelion
M: mean anomaly
''' 
def oe(r, r_dot, JD, a_jpl, e_jpl, i_jpl, O_jpl, w_jpl, M_jpl):

# iterates through position and velocity vectors, converts values into AU and AU/Gaussin Day, 
# and calculates position and velocity magnitudes
    for i in range(len(r)):
        r[i] = r[i]*(6.6845871222684*(10**-9))
    r_mag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)

    for i in range(len(r_dot)):
        r_dot[i] = r_dot[i]*(6.6845871222684*(10**-9))*(86400)*(58.13244087)
    r_dot_mag = np.sqrt(r_dot[0]**2 + r_dot[1]**2 + r_dot[2]**2)

# calculates orbital elements

# 1. a
    a = ((2/r_mag) - np.dot(r_dot, r_dot))**-1

# calculates angular momentum and angular momentum magnitude
    h = np.cross(r, r_dot)
    h_mag = np.sqrt(h[0]**2 + h[1]**2 + h[2]**2)

# 2. e
    e = np.sqrt(1 - ((h_mag**2)/a))

# 3. i
    i = np.arccos(h[2]/h_mag)

# 4. O

# first calculates arcsin and arccos values of O
    O_arcsin = np.arcsin(h[0]/(h_mag*np.sin(i)))
    O_arccos = np.arccos(-h[1]/(h_mag*np.sin(i)))

# performs quadrant check to determine true value of O
    if O_arcsin > 0:
        O = O_arccos
    else:
        O = -1*O_arccos
    O = O%(2*np.pi)

# 5. w

# first calculates U (distance from ascending node to asteroid) and v (true anomaly)

# calculates arcsin and arccos values for U
    U_arcsin = np.arcsin(r[2]/(r_mag*np.sin(i)))
    U_arccos = np.arccos((r[0]*np.cos(O) + r[1]*np.sin(O))/r_mag)

# performs quadrant check to determine true value of U
    if U_arcsin > O:
        U = U_arccos
    else:
        U = -1*U_arccos
    U = U%(2*np.pi)

# calculates arcsin and arccos values for v
    v_arcsin = np.arcsin = np.arcsin(((a*(1 - e**2))/(h_mag*e))*((np.dot(r, r_dot))/r_mag))
    v_arccos = np.arccos((a*(1-e**2) - r_mag)/(r_mag*e))

# performs quadrant check to determine true value of v
    if v_arcsin > 0:
        v = v_arccos
    else:
        v = -1*v_arccos
    v = v%(2*np.pi)

# calculates w
    w = (U - v)%(2*np.pi)

# 6. M

# first calculates E (eccentric anomaly)

# calculates arccos of E
    E_arccos = np.arccos((1/e)*(1 - (r_mag/a)))

# performs quadrant check to determine true value of E
    if np.sin(v) > 0:
        E = E_arccos
    else:
        E = -1*E_arccos
    E = E%(2*np.pi)

# calculates M
    M = (E - e*np.sin(E))%(2*np.pi)

# performs quadrant check to determine true value of M
    if np.sin(v) > np.pi and M < np.pi:
        M = (2*np.pi) - M
    elif v < np.pi and M > np.pi:
        M = (2*np.pi) - M
    else:
        M = M

# converts semimajor axis from JPL Horizons into AU
    a_jpl = a_jpl*(6.6845871222684*(10**-9))

# series of if statements to test the percent difference between calculated values and JPL Horizons values of orbital elements
    if ((abs(a - a_jpl)/((a + a_jpl)/2))*100):
        return ("semimajor axis inaccuracy")
    if ((abs(e - e_jpl)/((e + e_jpl)/2))*100):
        return ("eccentricity inaccuracy")
    if ((abs(i - i_jpl)/((i + i_jpl)/2))*100):
        return ("inclination inaccuracy")
    if ((abs(O - O_jpl)/((O + O_jpl)/2))*100):
        return ("right ascension of the ascending node inaccuracy")
    if ((abs(w - w_jpl)/((w + w_jpl)/2))*100):
        return ("argument of the perihelion inaccuracy")
    if ((abs(M - M_jpl)/((M + M_jpl)/2))*100):
        return ("mean anomaly inaccuracy")
    
    # given accuracy, returns a list of orbital elements
    else:
        return [a, e, rad_to_deg(i), rad_to_deg(O), rad_to_deg(w), M]


'''
calculates orbital elements using position and velocity vectors (same as "oe" function), but designed specifically for MoG
assumptions (all achieved in MoG):
1. position, velocity, and JPL Horizons orbital elements data have been calculated
2. position and velocity vectors have been calculated in AU and AU/Gaussian Day, respectively
'''
def oe_MoG(r, r_dot, a_jpl, e_jpl, i_jpl, O_jpl, w_jpl, M_jpl):

# calculates position and velocity magnitudes
    r_mag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
    r_dot_mag = np.sqrt(r_dot[0]**2 + r_dot[1]**2 + r_dot[2]**2)

# calculates orbital elements

# 1. a
    a = ((2/r_mag) - np.dot(r_dot, r_dot))**-1

# calculates angular momentum and angular momentum magnitude
    h = np.cross(r, r_dot)
    h_mag = np.sqrt(h[0]**2 + h[1]**2 + h[2]**2)

# 2. e
    e = np.sqrt(1 - ((h_mag**2)/a))

# 3. i
    i = np.arccos(h[2]/h_mag)

# 4. O

# first calculates arcsin and arccos values of O
    O_arcsin = np.arcsin(h[0]/(h_mag*np.sin(i)))
    O_arccos = np.arccos(-h[1]/(h_mag*np.sin(i)))

# performs quadrant check to determine true value of O
    if O_arcsin > 0:
        O = O_arccos
    else:
        O = -1*O_arccos
    O = O%(2*np.pi)

# 5. w

# first calculates U (distance from ascending node to asteroid) and v (true anomaly)

# calculates arcsin and arccos values for U
    U_arcsin = np.arcsin(r[2]/(r_mag*np.sin(i)))
    U_arccos = np.arccos((r[0]*np.cos(O) + r[1]*np.sin(O))/r_mag)

# performs quadrant check to determine true value of U
    if U_arcsin > O:
        U = U_arccos
    else:
        U = -1*U_arccos
    U = U%(2*np.pi)

# calculates arcsin and arccos values for v
    v_arcsin = np.arcsin = np.arcsin(((a*(1 - e**2))/(h_mag*e))*((np.dot(r, r_dot))/r_mag))
    v_arccos = np.arccos((a*(1-e**2) - r_mag)/(r_mag*e))

# performs quadrant check to determine true value of v
    if v_arcsin > 0:
        v = v_arccos
    else:
        v = -1*v_arccos
    v = v%(2*np.pi)

# calculates w
    w = (U - v)%(2*np.pi)

# 6. M

# first calculates E (eccentric anomaly)

# calculates arccos of E
    E_arccos = np.arccos((1/e)*(1 - (r_mag/a)))

# performs quadrant check to determine true value of E
    if np.sin(v) > 0:
        E = E_arccos
    else:
        E = -1*E_arccos
    E = E%(2*np.pi)

# calculates M
    M = (E - e*np.sin(E))%(2*np.pi)

# performs quadrant check to determine true value of M
    if np.sin(v) > np.pi and M < np.pi:
        M = (2*np.pi) - M
    elif v < np.pi and M > np.pi:
        M = (2*np.pi) - M
    else:
        M = M

# calculates percent difference between each orbital element and its corresponding JPL Horizons value
    a_diff = (abs(a - a_jpl)/((a + a_jpl/2)))*100
    e_diff = (abs(e - e_jpl)/((e + e_jpl/2)))*100
    i_diff = (abs(i - deg_to_rad(i_jpl))/((i + deg_to_rad(i_jpl/2))))*100
    O_diff = (abs(O - deg_to_rad(O_jpl))/((O + deg_to_rad(O_jpl/2))))*100
    w_diff = (abs(w - deg_to_rad(w_jpl))/((w + deg_to_rad(w_jpl/2))))*100
    M_diff = (abs(M - deg_to_rad(M_jpl))/((M + deg_to_rad(M_jpl/2))))*100

# prints orbital elements
    print("Orbital Elements:")
    print("Semimajor Axis: " + str(a))
    print("Eccentricity: " + str(e))
    print("Inclination: " + str(rad_to_deg(i)))
    print("Right Ascension of the Ascending Node: " + str(rad_to_deg(O)))
    print("Argument of the Perihelion: " + str(rad_to_deg(w)))
    print("Mean Anomaly: " + str(rad_to_deg(M)))

# prints percent differences
    print("Percent Differences between Calculated and JPL Horizons Values:")
    print("Semimajor Axis Percent Difference: " + str(a_diff))
    print("Eccentricity Percent Difference: " + str(e_diff))
    print("Inclination Percent Difference: " + str(i_diff))
    print("Right Ascension of the Ascending Node Percent Difference: " + str(O_diff))
    print("Argument of the Perihelion Percent Difference: " + str(w_diff))
    print("Mean Anomaly Percent Difference: " + str(M_diff))

# returns a list of orbital elements
    return [a, e, rad_to_deg(i), rad_to_deg(O), rad_to_deg(w), M]


# given a text file appropriately formatted for the "ephemeris_gen" function, extracts data
def ephemeris_gen_data(file):

# opens text file and creates a list of lines
    data = open(file)
    content = data.readlines()

# extracts each orbital element from the list, converts each from string to float, and converts units as follows:
# a: km to au
# i, O, w, M: degrees to radians
    a = float(content[0])/149597871
    e = float(content[1])
    i = deg_to_rad(float(content[2]))
    O = deg_to_rad(float(content[3]))
    w = deg_to_rad(float(content[4]))
    M_ref = deg_to_rad(float(content[5]))

# extracts reference and target JD values and converts them into floats
    ref_JD = float(content[6])
    t_JD = float(content[7])

# extracts reference RA/dec values in hms/dms, respectively, as strings
# splits strings, converts strings into floats, and converts hms/dms into degrees
    ref_RA = content[8]
    ref_dec = content[9]
    ref_RA_split = ref_RA.split(' ')
    ref_dec_split = ref_dec.split(' ')
    ref_RA_deg = (float(ref_RA_split[0]) + (float(RA[1])/60) + (float(ref_RA_split[2])/3600))*15
    ref_dec_deg = (float(ref_dec_split[0]) + (float(ref_dec_split[1])/60) + (float(ref_dec_split[2])/3600))

# extracts the Earth/Sun vector as a string, splits it into its components, and organizes it in an array
    es_str = content[10]
    es_split = es_str.split(' ')
    es = np.array([[float(es_split[0])], [float(es_split[1])], [float(es_split[2])]])

# returns orbital elements, JD values, RA/dec values, and the Earth/Sun vector
    return [a, e, i, O, w, M_ref, ref_JD, t_JD, ref_RA, ref_dec, es]


'''
uses orbital elements at a specific date, reference and target times, and the Earth/Sun vector
to calculate RA/dec values and their respective uncertainties
assumptions (all achieved by "ephemeris_gen_data" function):
1. orbital elements are provided with proper units (a: AU, i, O, w, M_ref: degrees)
2. RA/dec are floats and provided in degrees
3. Earth/Sun vector is a vector stored as an np.array
'''
def ephemeris_gen(a, e, i, O, w, M_ref, ref_JD, t_JD, ref_RA, ref_dec, es):
    
# calculates orbital period and calculates Gaussian Days into days
    T = np.sqrt(4*(np.pi**2)*(a**3))*58.13244087

# calculates JD of the time at which the object last passed the perihelion distance
    t_peri = -((M_ref*(T)/(2*np.pi)) - ref_JD)

# calculates mean anomaly at the target JD
    M = ((2*np.pi*(t_JD - t_peri))/(T))%(2*np.pi)

# performs Newton's Method to calculate the eccentric anomaly
    E0 = M
    for n in range(100000):
        E0 = E0 - ((M - (E0 - e*np.sin(E0)))/(-1 + e*np.cos(E0)))
    E = E0%(2*np.pi)

# calculates the orbital position components, and stores the coordinates in a 1x3 array
    x_orb = a*np.cos(E) - a*e
    y_orb = a*(np.sqrt(1 - (e**2)))*np.sin(E)
    z_orb = 0
    orb_coords = np.array([[x_orb], [y_orb], [z_orb]])

# stores rotation matrices (used to convert to ecliptic coordinates) in arrays
    w_rot = np.array([[np.cos(w), -np.sin(w), 0], [np.sin(w), np.cos(w), 0], [0, 0, 1]])
    i_rot = np.array([[1, 0, 0], [0, np.cos(i), -np.sin(i)], [0, np.sin(i), np.cos(i)]])
    O_rot = np.array([[np.cos(O), -np.sin(O), 0], [np.sin(O), np.cos(O), 0], [0, 0, 1]])

# performs matrix multiplication to convert orbital coordinates to ecliptic coordinates
    ecl_coords = O_rot@i_rot@w_rot@orb_coords

# creates a rotation matrix to convert ecliptic coordinates to equatorial coordinates
    tilt = deg_to_rad(23.44)
    ecl_eq_rot = np.array([1, 0, 0], [0, np.cos(tilt), -np.sin(tilt)], [0, np.sin(tilt), np.cos(tilt)])

# performs matrix multiplication to convert equatorial coordinates to ecliptic coordinates
    eq_coords = ecl_eq_rot@ecl_coords

# first converts equatorial coordinates from AU to km, then calculates rho (Earth/object vector)
# by adding Earth/Sun vector and equatorial coordinates
    rho = es + (eq_coords*149597871)

# extracts rho components
    rho_x = rho[0]
    rho_y = rho[1]
    rho_z = rho[2]

# calculates rho magnitude
    rho_mag = np.sqrt((rho_x**2) + (rho_y**2) + (rho_z**2))

# calculates RA/dec and converts them to degrees
    RA = rad_to_deg(float(np.atan2(rho_y, rho_x)%(2*np.pi)))
    dec = rad_to_deg(float(np.arcsin(rho_z/rho_mag)%(2*np.pi)))

# calculates percent difference (in seconds) between calculated RA/dec and reference RA/dec
    RA_pd = 3600*(abs(RA - ref_RA))/(ref_RA)*100
    dec_pd = 3600*(abs(dec - dec))/(ref_dec)

# prints reference and calculated values of RA/dec and the percent difference for each
    print("Reference RA/dec: " + str(ref_RA) + "," + str(ref_dec))
    print("RA/dec: " + str(RA) + "," + str(dec))
    print("Percent difference in RA/dec: " + str(RA_pd) + "/" + str(dec_pd))

# returns list of RA/dec
    return [RA, dec]


'''
calculates the f and g constants for MoG using a Taylor approximation:
r, r_dot are position, velocity vectors
tau is the time interval
order4 is a boolean value that determines whether the Taylor approximation will be 3rd or 4th order
'''
def fg(r, r_dot, tau, order4):
    
# calculates position and velocity magnitudes
    r_mag = (float(r[0])**2 + float(r[1])**2 + float(r[2])**2)**0.5
    v_mag = (float(r_dot[0])**2 + float(r_dot[1])**2 + float(r_dot[2])**2)**0.5

# calculates u, z, q values
    u = 1/(r_mag**3)
    z = float(np.dot(r, r_dot)/(r_mag**2))
    q = v_mag**2/(r_mag**2) - u

# based on boolean value of order4, performs 3rd or 4th order Taylor approximation
    if order4 == True:
        f = 1 - (1/2)*u*(float(tau)**2) + (1/2)*(u*z)*(float(tau)**2) + (1/24)*(3*u*q - 15*u*(z**2) + (u**2))*(float(tau)**4)
        g = (float(tau)) - (1/6)*u*(float(tau)**3) + (1/4)*(u*z)*(float(tau)**4)
    else:
        f = 1 - (1/2)*u*(float(tau)**2) + (1/2)*(u*z)*(float(tau)**2)
        g = (float(tau)) - (1/6)*u*(float(tau)**3)

# returns a list of the f and g constants
    return ([f, g])


'''
The following function performs the Method of Gauss (MoG)

Inputs:
1. times of 3 observations (provided in JD)
2. RA/dec of the object at those 3 times (provided in radians)
3. Earth/Sun vector at those 3 times (provided in AU)
4. JPL Horizons reference orbital elements

Outputs:
1. orbital elements
'''
def MoG(t1, t2, t3, RA1, RA2, RA3, dec1, dec2, dec3, es1, es2, es3, a_ref, e_ref, i_ref, O_ref, w_ref, M_ref):

# SETUP

# declares k (constant) and initial tau (time interval) values
    k = 0.0172020989484
    tau1 = k*(t1 - t2)
    tau3 = k (t3 - t2)
    tau0 = k*(t3 - t1)

# FIRST ITERATION (based on using tau values to "guess" initial values)

# calculates c-constants
    c1 = tau3/tau0
    c2 = -1
    c3 = -tau1/tau0

# calculates rho unit vectors
    rho1_hat = np.array([float(np.cos(RA1)*np.cos(dec1)),float(np.sin(RA1)*np.cos(dec1)),float(np.sin(dec1))])
    rho2_hat = np.array([float(np.cos(RA2)*np.cos(dec2)),float(np.sin(RA2)*np.cos(dec2)),float(np.sin(dec2))])
    rho3_hat = np.array([float(np.cos(RA3)*np.cos(dec3)),float(np.sin(RA3)*np.cos(dec3)),float(np.sin(dec3))])
    

# calculates D-constants

    D0 = float(rho1_hat@(np.cross(rho2_hat,rho3_hat)))

    D11 = float((np.cross(es1,rho2_hat))@rho3_hat)
    D12 = float((np.cross(es2,rho2_hat))@rho3_hat)
    D13 = float((np.cross(es3,rho2_hat))@rho3_hat)

    D21 = float((np.cross(rho1_hat,es1))@rho3_hat)
    D22 = float((np.cross(rho1_hat,es2))@rho3_hat)
    D23 = float((np.cross(rho1_hat,es3))@rho3_hat)

    D31 = float(rho1_hat@(np.cross(rho2_hat,es1)))
    D32 = float(rho1_hat@(np.cross(rho2_hat,es2)))
    D33 = float(rho1_hat@(np.cross(rho2_hat,es3)))
    
# uses c-constants and D-constants to calculate initial rho-magnitudes
    rho1_mag = (c1*D11+c2*D12+c3*D13)/(c1*D0)
    rho2_mag = (c1*D21+c2*D22+c3*D23)/(c2*D0)
    rho3_mag = (c1*D31+c2*D32+c3*D33)/(c3*D0)

# declares speed of light in AU/day and sets while-loop counter to 0
    c =  173.144643267
    count = 0

# performs first light travel time correction
    t1_corr = t1-(rho1_mag/c)
    t2_corr = t2-(rho2_mag/c)
    t3_corr = t3-(rho3_mag/c)

# recalculates tau values using corrected time values
    tau1 = k*(t1_corr-t2_corr)
    tau3 = k*(t3_corr-t2_corr)
    tau0 = k*(t3_corr-t1_corr)

# performs first calculation of position vectors (in eq plane) using fundamental triangle
    r1 = rho1_mag*rho1_hat - es1
    r2 = rho2_mag*rho2_hat - es2
    r3 = rho3_mag*rho3_hat - es3

# performs first calculation of average velocities between observation times
    r12_dot = (r2-r1)/(-tau1)
    r23_dot = (r3-r2)/(tau3)

# performs first calculation of velocity vector using a weighted average of average velocities
    r2_dot = ((tau3/tau0)*r12_dot) - ((tau1/tau0)*r23_dot)

# SUBSEQUENT ITERATIONS (while-loop based on Taylor approximation to accurately calculate position/velocity vectors)

# asserts check as r2 (will be optimized during the while-loop)
# and change as 10 (any n>10**-12 will suffice as it will be great enough to initiate the while-loop)
    check = r2
    change = 10

# calculates initial f- and g-constants using "fg" function
    f1 = fg(r2,r2_dot,tau1,True)[0]
    g1 = fg(r2,r2_dot,tau1,True)[1]
    f3 = fg(r2,r2_dot,tau3,True)[0]
    g3 = fg(r2,r2_dot,tau3,True)[1]

# while loop runs until the r2 value converges or the count variable exceeds 1000
    while change>10**-12 and count<1000:

# recalculates corrected time values using recalculated rho-magnitude values 
        t1_corr = t1 - (rho1_mag/c)
        t2_corr = t2 - (rho2_mag/c)
        t3_corr = t3 - (rho3_mag/c)

# recalculates tau values using recalculated time values
        tau1 = k*(t1_corr-t2_corr)
        tau3 = k*(t3_corr-t2_corr)
        tau0 = k*(t3_corr-t1_corr)

# recalculates f- and g-constants using recalculated tau values
        f1 = fg(r2,r2_dot,tau1,True)[0]
        g1 = fg(r2,r2_dot,tau1,True)[1]
        f3 = fg(r2,r2_dot,tau3,True)[0]
        g3 = fg(r2,r2_dot,tau3,True)[1]

# recalculates c- and d-constants using recalculated f- and g-constants
        c1 = (g3)/(f1*g3-g1*f3)
        c3 = (-g1)/(f1*g3-g1*f3)
        d1 = (-f3)/(f1*g3-f3*g1)
        d3 = (f1)/(f1*g3-f3*g1)

# recalculates rho-magnitudes using D-constants and recalculated c-constants
        rho1_mag = (c1*D11+c2*D12+c3*D13)/(c1*D0)
        rho2_mag = (c1*D21+c2*D22+c3*D23)/(c2*D0)
        rho3_mag = (c1*D31+c2*D32+c3*D33)/(c3*D0)

# recalculates position vectors using rho unit vectors, Earth/Sun vectors, and recalculated rho-magnitudes
        r1 = rho1_mag*rho1_hat - es1
        r2 = rho2_mag*rho2_hat - es2
        r3 = rho3_mag*rho3_hat - es3

# recalculaes velocity vector using d-constants and recalculated position vectors
        r2_dot = d1*r1 + d3*r3


# recalculates change using recalculated r2 position vector
# asserts check as the recalculated r2 position vector
# adds 1 to the count
        change = np.linalg.norm(r2-check)/np.linalg.norm(r2)
        check = r2
        count += 1
    
# rotates position/velocity vectors from equatorial to ecliptic plane by asserting the tilt angle, creating the rotation matrix, 
# and applying it to the position and velocity vectors
    tilt = deg_to_rad(23.44)
    eq_to_ecl = np.array([[1, 0, 0], [0, np.cos(tilt), np.sin(tilt)], [0, -np.sin(tilt), np.cos(tilt)]])
    r2_ecl = eq_to_ecl@r2
    r2_dot_ecl = eq_to_ecl@r2_dot

# 1. uses the "oe_MoG" function and the position/velocity vectors to calculate the orbital elements and percent difference values
# between the calculated orbital elements and the JPL Horizons values
# 2. prints the orbital elements and the percent difference values
# 3. returns a list of the orbital elements
    return oe_MoG(r2_ecl, r2_dot_ecl, a_ref, e_ref, i_ref, O_ref, w_ref, M_ref)

# calculates the uncertainty associated with each calculated orbital element value
def get_uncertainty(t1, t2, t3, RA1, RA1_uncertainty, RA2, RA2_uncertainty, RA3, RA3_uncertainty, dec1, dec1_uncertainty, dec2, dec2_uncertainty, dec3, dec3_uncertainty, N, es1, es2, es3, a_ref, e_ref, i_ref, O_ref, w_ref, M_ref):
    
# creates arrays to append orbital element values to
    a_values = np.array([])
    e_values = np.array([])
    i_values = np.array([])
    O_values = np.array([])
    w_values = np.array([])
    M_values = np.array([])

# for-loop runs 'N' times
# calculates orbital elements using RA/dec values randomly selected from their Gaussian distributions
    for i in range(N):
       
# uses the uncertainty values of RA/dec to select random RA/dec values from their Gaussian distributions
        RA1_u = np.random.normal(RA1, RA1_uncertainty)
        RA2_u = np.random.normal(RA2, RA2_uncertainty)
        RA3_u = np.random.normal(RA3, RA3_uncertainty)
        dec1_u = np.random.normal(dec1, dec1_uncertainty)
        dec2_u = np.random.normal(dec2, dec2_uncertainty)
        dec3_u = np.random.normal(dec3, dec3_uncertainty)

# uses the "MoG" function and random RA/dec values to calculate orbital elements
        orb_elements = MoG(t1, t2, t3, RA1_u, RA2_u, RA3_u, dec1_u, dec2_u, dec3_u, es1, es2, es3, a_ref, e_ref, i_ref, O_ref, w_ref, M_ref)
        
# extracts orbital elements from the list returned by "MoG" function
        a = orb_elements[0]
        e = orb_elements[1]
        i = orb_elements[2]
        O = orb_elements[3]
        w = orb_elements[4]
        M = orb_elements[5]

# appends orbital elements to their respective arrays
        a_values = np.append(a_values, a)
        e_values = np.append(e_values, e)
        i_values = np.append(i_values, i)
        O_values = np.append(O_values, O)
        w_values = np.append(w_values, w)
        M_values = np.append(M_values, M)

# calculates uncertainty of each orbital element by calculating the standard deviation of each list of orbital elements
    a_uncertainty = np.std(a_values)
    e_uncertainty = np.std(e_values)
    i_uncertainty = np.std(i_values)
    O_uncertainty = np.std(O_values)
    w_uncertainty = np.std(w_values)
    M_uncertainty = np.std(M_values)

# prints uncertainty values
    print("Uncertainty Values:")
    print(a_uncertainty)
    print(e_uncertainty)
    print(i_uncertainty)
    print(O_uncertainty)
    print(w_uncertainty)
    print(M_uncertainty)       
    
# returns list of uncertainty values
    return [a_uncertainty, e_uncertainty, i_uncertainty, O_uncertainty, w_uncertainty, M_uncertainty]


# given RA/dec, RA/dec uncertainty, JPL Horizons RA/dec reference values, and an integer 'N' to determine for-loop iterations,
# plots the Gaussian distribution of RA/dec
def RA_dec_hist(RA, RA_u, RA_ref, dec, dec_u, dec_ref, N):

# prepares empty arrays to append RA/dec values to
    RA_values = np.array([])
    dec_values = np.array([])

# for-loop runs 'N' times
# selects random RA/dec using uncertainty values and appends them to arrays
    for i in range(N):
        RA_append = np.random.normal(RA, RA_u)
        RA_values = np.append(RA_values, RA_append)
        dec_append = np.random.normal(dec, dec_u)
        dec_values = np.append(dec_values, dec_append)
    
# plots RA/dec histograms with:
# 1. vertical lines indicating observed values, uncertainty bounds, and reference values
# 2. x and y axis labels
# 3. legends
# 4. plot titles

# RA
    plt.hist(RA_values)
    plt.axvline(x=RA, color='purple', label='observed value', linewidth=5)
    plt.axvline(x=RA+RA_u, color='red', label='upper uncertainty bound')
    plt.axvline(x=RA-RA_u, color='red', label='lower uncertainty bound')
    plt.axvline(x=RA_ref, color='pink', label='JPL reference value', linewidth=5)
    plt.xlabel('value')
    plt.ylabel('count')
    plt.legend(loc='upper right')
    plt.title('RA')
    plt.show()

# dec
    plt.hist(dec_values)
    plt.axvline(x=dec,color='purple',label='observed value', linewidth=5)
    plt.axvline(x=dec+dec_u,color='red',label='upper uncertainty bound')
    plt.axvline(x=dec-dec_u,color='red',label='lower uncertainty bound')
    plt.axvline(x=dec_ref,color='pink',label='JPL reference value',linewidth=5)
    plt.xlabel('value')
    plt.ylabel('count')
    plt.legend(loc='upper right')
    plt.title('Declination')
    plt.show()

    return

# given semimajor axis, semimajor axis uncertainty value, JPL Horizons semimajor axis reference value, 
# and an integer 'N' to determine for-loop iterations, plots the Gaussian distribution of semimajor axis 
def a_hist(a, a_u, a_ref, N):

# prepares empty array to append a values to
    a_values = np.array([])

# for-loop runs 'N' times
# selects random a values using uncertainty value and appends them to the array
    for n in range(N):
        a_append = np.random.normal(a, a_u)
        a_values = np.append(a_values, a_append)
    
# plots a histogram with:
# 1. vertical lines indicating calculated value, uncertainty bounds, and reference value
# 2. x and y axis labels
# 3. a legend
# 4. a plot title
    plt.hist(a_values)
    plt.axvline(x=a, color='purple', label='calculated value', linewidth=5)
    plt.axvline(x=a+a_u, color='red', label='upper uncertainty bound')
    plt.axvline(x=a-a_u, color='red', label='lower uncertainty bound')
    plt.axvline(x=a_ref, color='pink', label='JPL reference value', linewidth=5)
    plt.xlabel('value')
    plt.ylabel('count')
    plt.legend(loc='upper right')
    plt.title('a')
    plt.show()  

    return

# given eccentricity, eccentricity uncertainty value, JPL Horizons eccentricity reference value, 
# and an integer 'N' to determine for-loop iterations, plots the Gaussian distribution of eccentricity 
def e_hist(e, e_u, e_ref, N):

# prepares empty array to append e values to
    e_values = np.array([])

# for-loop runs 'N' times
# selects random a values using uncertainty value and appends them to the array
    for n in range(N):
        e_append = np.random.normal(e, e_u)
        e_values = np.append(e_values, e_append)
    
# plots a histogram with:
# 1. vertical lines indicating calculated value, uncertainty bounds, and reference value
# 2. x and y axis labels
# 3. a legend
# 4. a plot title
    plt.hist(e_values)
    plt.axvline(x=e, color='purple', label='calculated value', linewidth=5)
    plt.axvline(x=e+e_u, color='red', label='upper uncertainty bound')
    plt.axvline(x=e-e_u, color='red', label='lower uncertainty bound')
    plt.axvline(x=e_ref, color='pink', label='JPL reference value', linewidth=5)
    plt.xlabel('value')
    plt.ylabel('count')
    plt.legend(loc='upper right')
    plt.title('e')
    plt.show()  

    return

# given inclination, inclination uncertainty value, JPL Horizons inclination reference value, 
# and an integer 'N' to determine for-loop iterations, plots the Gaussian distribution of inclination 
def i_hist(i, i_u, i_ref, N):

# prepares empty array to append e values to
    i_values = np.array([])

# for-loop runs 'N' times
# selects random a values using uncertainty value and appends them to the array
    for n in range(N):
        i_append = np.random.normal(i, i_u)
        i_values = np.append(i_values, i_append)
    
# plots a histogram with:
# 1. vertical lines indicating calculated value, uncertainty bounds, and reference value
# 2. x and y axis labels
# 3. a legend
# 4. a plot title
    plt.hist(i_values)
    plt.axvline(x=i, color='purple', label='calculated value', linewidth=5)
    plt.axvline(x=i+i_u, color='red', label='upper uncertainty bound')
    plt.axvline(x=i-i_u, color='red', label='lower uncertainty bound')
    plt.axvline(x=i_ref, color='pink', label='JPL reference value', linewidth=5)
    plt.xlabel('value')
    plt.ylabel('count')
    plt.legend(loc='upper right')
    plt.title('e')
    plt.show()  

    return

# given right ascension of the ascending node, right ascension of the ascending node uncertainty value, 
# JPL Horizons right ascension of the ascending node reference value, and an integer 'N' to determine for-loop iterations, 
# plots the Gaussian distribution of right ascension of the ascending node 
def O_hist(O, O_u, O_ref, N):

# prepares empty array to append e values to
    O_values = np.array([])

# for-loop runs 'N' times
# selects random a values using uncertainty value and appends them to the array
    for n in range(N):
        O_append = np.random.normal(O, O_u)
        O_values = np.append(O_values, O_append)
    
# plots a histogram with:
# 1. vertical lines indicating calculated value, uncertainty bounds, and reference value
# 2. x and y axis labels
# 3. a legend
# 4. a plot title
    plt.hist(O_values)
    plt.axvline(x=O, color='purple', label='calculated value', linewidth=5)
    plt.axvline(x=O+O_u, color='red', label='upper uncertainty bound')
    plt.axvline(x=O-O_u, color='red', label='lower uncertainty bound')
    plt.axvline(x=O_ref, color='pink', label='JPL reference value', linewidth=5)
    plt.xlabel('value')
    plt.ylabel('count')
    plt.legend(loc='upper right')
    plt.title('O')
    plt.show()  

    return

# given argument of the perihelioin, argument of the perihelion uncertainty value, 
# JPL Horizons argument of the perihelion reference value, # and an integer 'N' to determine for-loop iterations, 
# plots the Gaussian distribution of argument of the perihelion 
def w_hist(w, w_u, w_ref, N):

# prepares empty array to append e values to
    w_values = np.array([])

# for-loop runs 'N' times
# selects random a values using uncertainty value and appends them to the array
    for n in range(N):
        w_append = np.random.normal(w, w_u)
        w_values = np.append(w_values, w_append)
    
# plots a histogram with:
# 1. vertical lines indicating calculated value, uncertainty bounds, and reference value
# 2. x and y axis labels
# 3. a legend
# 4. a plot title
    plt.hist(w_values)
    plt.axvline(x=w, color='purple', label='calculated value', linewidth=5)
    plt.axvline(x=w+w_u, color='red', label='upper uncertainty bound')
    plt.axvline(x=w-w_u, color='red', label='lower uncertainty bound')
    plt.axvline(x=w_ref, color='pink', label='JPL reference value', linewidth=5)
    plt.xlabel('value')
    plt.ylabel('count')
    plt.legend(loc='upper right')
    plt.title('w')
    plt.show()  

    return

# given mean anomaly, mean anomaly uncertainty value, JPL Horizons mean anomaly reference value, 
# and an integer 'N' to determine for-loop iterations, plots the Gaussian distribution of mean anomaly 
def M_hist(M, M_u, M_ref, N):

# prepares empty array to append e values to
    M_values = np.array([])

# for-loop runs 'N' times
# selects random a values using uncertainty value and appends them to the array
    for n in range(N):
        M_append = np.random.normal(M, M_u)
        M_values = np.append(M_values, M_append)
    
# plots a histogram with:
# 1. vertical lines indicating calculated value, uncertainty bounds, and reference value
# 2. x and y axis labels
# 3. a legend
# 4. a plot title
    plt.hist(M_values)
    plt.axvline(x=M, color='purple', label='calculated value', linewidth=5)
    plt.axvline(x=M+M_u, color='red', label='upper uncertainty bound')
    plt.axvline(x=M-M_u, color='red', label='lower uncertainty bound')
    plt.axvline(x=M_ref, color='pink', label='JPL reference value', linewidth=5)
    plt.xlabel('value')
    plt.ylabel('count')
    plt.legend(loc='upper right')
    plt.title('M')
    plt.show()  

    return

