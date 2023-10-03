import re
import csv
import math
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord


def main():
    g_data, r_data = file_controller()
    key = read_key()
    # add v-magnitude to key split file
    # ref = read_file_ref()
    # vm = read_file_vm(ref)
    # key_vmag = map_key_mag(key, vm)
    # histogram(key_vmag)
    # csv_key_export(key)

    sample_g, sample_r, center_point = sample_collecting(g_data, r_data)
    # print("SAMPLE G IS ", len(sample_g), "\n", sample_g )
    # print("SAMPLE R IS ", len(sample_r), "\n", sample_r)
    # print("CENTER POINTS IS ", len(center_point), "\n",  center_point)
    all_delta_ra_dec = difference(sample_g, sample_r)
    data_shifting(r_data, all_delta_ra_dec)
    # # time to rotate
    G = data_rotate(g_data, r_data, center_point, -0.40)

    # # theta_controller(G, R, center_point)

    matched = matching(key, G)
    plot_key_data(G, key, 'The plot shows all G data and key')
    # print(matched)
    plot_filnal(matched, 'The plot shows matched data and key')
    csv_export(matched)


def file_controller():
    data_g, data_r = [], []
    header = "dataset\\"
    g_file = ['g23_n.cat', 'g32_n.cat']
    r_file = ['r23_n.cat', 'r32_n.cat']

    for i in range(len(g_file)):
        g_list = read_file(header + g_file[i])
        r_list = read_file(header + r_file[i])
        # filter magnitude(g&r) in range 15 - 25 of filter G & R
        g_updated = filter_magnitude(g_list)
        r_updated = filter_magnitude(r_list)
        data_g.append(g_updated)
        data_r.append(r_updated)
    return data_g, data_r


def read_file(filename):
    data = []
    myfile = open(filename)
    info = myfile.readlines()
    for i in range(len(info)):
        info[i].rstrip('\n')
        a = re.findall(r"\S+", info[i])
        if(a[0] != '#'):
            temp = []
            temp.append(float(a[12]))   # 0: RA
            temp.append(float(a[13]))   # 1: DEC
            temp.append(a[1])           # 2: MAG_APER
            temp.append(filename)       # 3: file name
            temp.append(a[0])           # 4: Object number
            temp.append(float(a[14]))   # 5: X
            temp.append(float(a[15]))   # 6: Y
            data.append(temp)
        myfile.close()
    return data


def read_key():
    filename = 'aj341087t3_ascii.txt'
    header = "dataset\\"
    myfile = open(header + filename)
    info = myfile.readlines()
    key_splited = []
    for j in range(len(info)):
        temp = []
        if(j >= 5):  # start at line 5 [you can change parameter here.]
            info[j] = info[j].rstrip('\n')
            info[j] = re.findall(r'\S+', info[j])
            if filename == 'aj341087t2_ascii.txt':
                temp.append(info[j][4])  # get RA
                temp.append(info[j][5])  # get Dec
            elif filename == 'aj341087t3_ascii.txt':
                ra = info[j][4] + ':' + info[j][5] + ':' + info[j][6]
                dec = info[j][7] + ':' + info[j][8] + ':' + info[j][9]
                temp.append(ra)  # get RA
                temp.append(dec)  # get DEC
            temp.append(info[j][0])  # get No.
            temp.append(info[j][2])  # get PR95
            key_splited.append(temp)
            myfile.close()
    # convert RA and DEC to floating form
    key = convert_to_degree(key_splited)
    return key

def map_key_mag(ra_dec_K, vm):
    key_vmag = []
    for j in range(len(ra_dec_K)):
        temp = []
        for k in range(len(vm)):
            # if reference number are matched with ID of key file
            if vm[k][1] == ra_dec_K[j][2] or vm[k][2] == ra_dec_K[j][3]:
                # then add v_mag into key file in index at ID
                temp.append(ra_dec_K[j][0])  # RA
                temp.append(ra_dec_K[j][1])  # DEC
                temp.append(vm[k][3])   # V_mag
        if temp != []:
            key_vmag.append(temp)
    return key_vmag


def read_file_ref():
    splited_ref = []
    ref = []
    myfile = open("dataset\\aj403657t2_mrt[3791].txt")
    info = myfile.readlines()
    for i in range(29, len(info)):  # start split data at line 29
        info[i] = info[i].rstrip('\n')
        info[i] = re.findall(r'\S+', info[i])
        splited_ref.append(info[i])
    myfile.close()
    for i in range(len(splited_ref)):
        temp = []
        temp.append(splited_ref[i][0])  # add number of line
        temp.append(splited_ref[i][7])  # add number of reference ID
        if len(splited_ref[i]) < 9:  # add PR95
            temp.append('0')
        else:
            temp.append(splited_ref[i][8])
        ref.append(temp)
    return ref


def read_file_vm(ref):
    splited_vm = []
    myfile = open('dataset\\aj403657t3_mrt[3788].txt')
    info = myfile.readlines()
    for i in range(32, len(info)):   # start split data at line 32
        info[i] = info[i].rstrip('\n')
        info[i] = re.findall(r'\S+', info[i])
        splited_vm.append(info[i])
    myfile.close()
    for i in range(len(splited_vm)):
        ref[i].append(float(splited_vm[i][1]))  # add v_mag
    return ref


def convert_to_degree(key_splited):
    for i in range(len(key_splited)):
        # use for file 'aj341087t2_ascii.txt' that RA and Dec have ':' between number
        ra = key_splited[i][0].split(':')
        dec = key_splited[i][1].split(':')
        # calculate HH:MM:SS to degree RA
        key_splited[i][0] = (float(ra[0]) + float(ra[1]) /
                             60 + float(ra[2])/3600)*15        # RA
        key_splited[i][1] = float(
            dec[0]) + float(dec[1])/60 + float(dec[2])/3600          # DEC

        # shift key condition
        key_splited[i][0] += 0.0113
        key_splited[i][1] -= 0.00276
        # key_splited[i][0] += 0.0
        # key_splited[i][1] -= 0.0
    return key_splited


def filter_magnitude(data):
    result = []
    for i in range(len(data)):
        # convert_mag_aper_to_gmagnitude
        data[i][2] = convert_to_mag(data[i][2])
        # if gmagnitude is in rangr (15-25)
        if data[i][2] >= 16 and data[i][2] <= 24:
            result.append(data[i])  # keep all data in set (RA,DEC,G_MAG)
    # sort the data by using magnitude. then select 30 first stars to the process
    sort_list = sorted(result, key=lambda l: l[2], reverse=False)
    # select top 30 row that are sorted
    result = sort_list[0:600]
    return result


def convert_to_mag(data):
    zeropoint = 26.5520
    exptime = 660.252
    result = float(data) + zeropoint + 2.5*(math.log(exptime, 10))
    return result


def sample_collecting(g_data, r_data):
    print(g_data)
    center_p_collect = []
    all_center_p_collect = []
    # collecting sample g&r by using PIXEL x = 1006-1106 and y = 2272-2372
    sample_g, sample_r = [], []
    for i in range(len(g_data)):
        center_point = []
        temp = []
        # min
        min_ra_g = min(g_data[i], key=lambda l: l[0])
        min_dec_g = min(g_data[i], key=lambda l: l[1])
        min_ra_r = min(r_data[i], key=lambda l: l[0])
        min_dec_r = min(r_data[i], key=lambda l: l[1])
        # max
        max_ra_g = max(g_data[i], key=lambda l: l[0])
        max_dec_g = max(g_data[i], key=lambda l: l[1])
        max_ra_r = max(r_data[i], key=lambda l: l[0])
        max_dec_r = max(r_data[i], key=lambda l: l[1])
        # set treshold
        ra_g_center = (max_ra_g[0] + min_ra_g[0])/2
        dec_g_center = (max_dec_g[1] + min_dec_g[1])/2
        ra_r_center = (max_ra_r[0] + min_ra_r[0])/2
        dec_r_center = (max_dec_r[1] + min_dec_r[1])/2
        temp.append(ra_g_center)
        temp.append(dec_g_center)
        temp.append(ra_r_center)
        temp.append(dec_r_center)
        all_center_p_collect.append(temp)
        # print(ra_g_center, dec_g_center, ra_r_center, dec_r_center)
        # center_point_of_image_using_g_center
        center_point.append(ra_g_center)
        center_point.append(dec_g_center)
        center_p_collect.append(center_point)
        print("Center Point is : ", center_point)

    for i in range(len(g_data)):
        temp_sample = []
        for j in range(len(g_data[i])):
            if g_data[i][j][0] >= all_center_p_collect[i][0]-0.05 and g_data[i][j][0] <= all_center_p_collect[i][0]+0.05 and g_data[i][j][1] >= all_center_p_collect[i][1]-0.05 and g_data[i][j][1] <= all_center_p_collect[i][1]+0.05:
                temp_sample.append(g_data[i][j])
        sample_g.append(temp_sample)
    for i in range(len(r_data)):
        temp_sample = []
        for j in range(len(r_data[i])):
            if r_data[i][j][0] >= all_center_p_collect[i][2]-0.05 and r_data[i][j][0] <= all_center_p_collect[i][2]+0.05 and r_data[i][j][1] >= all_center_p_collect[i][3]-0.05 and r_data[i][j][1] <= all_center_p_collect[i][3]+0.05:
                temp_sample.append(r_data[i][j])
        sample_r.append(temp_sample)
    return sample_g, sample_r, center_p_collect


def difference(sample_g, sample_r):
    number = [22, 23, 24, 31, 32, 33]
    all_delta_ra_dec = []
    for i in range(len(sample_g)):
        text = "The plot shows sample of coordiante in G & R " + str(number[i])
        temp = []
        plot_sample(sample_g[i], sample_r[i], text)
        print("***** Enter the Coordinate (RA0 and DEC0) of Sample G",
              number[i], "*****")
        ra0, dec0 = input().split()
        print("***** Enter the Coordinate (RA0 and DEC0) of Sample R",
              number[i], "*****")
        ra, dec = input().split()
        ra, ra0, dec, dec0 = float(ra), float(ra0), float(dec), float(dec0)
        # find the delta of each RA, DEC
        print(ra, ra0, dec, dec0)
        delta_ra = ra - ra0
        delta_dec = dec - dec0
        print("Delta_RA: ", delta_ra)
        print("Delta_DEC: ", delta_dec)
        # convert to shift constant
        delta_ra = delta_ra*(-1)
        delta_dec = delta_dec*(-1)
        temp.append(delta_ra)
        temp.append(delta_dec)
        print("Delta_RA: ", delta_ra)
        print("Delta_DEC: ", delta_dec)
        all_delta_ra_dec.append(temp)

        new_list_r = [[[0 for i in range(7)] for j in range(
            len(sample_r[i]))] for k in range(len(sample_r))]
        for j in range(len(sample_r[i])):
            # plus shift constant for coordinate ra of sample r
            new_list_r[i][j][0] = new_list_r[i][j][0] + \
                sample_r[i][j][0] + delta_ra
            # plus shift constant for coordinate dec of sample r
            new_list_r[i][j][1] = new_list_r[i][j][1] + \
                sample_r[i][j][1] + delta_dec
            new_list_r[i][j][2] = sample_r[i][j][2]
            new_list_r[i][j][3] = sample_r[i][j][3]
            new_list_r[i][j][4] = sample_r[i][j][4]
            new_list_r[i][j][5] += sample_r[i][j][5]
            new_list_r[i][j][6] += sample_r[i][j][6]
        # plot check
        text2 = "The plot shows shifted sample in G & R " + str(number[i])
        plot_sample(sample_g[i], new_list_r[i], text2)

        count = 0
        G, R = [], []
        for j in range(len(sample_g[i])):
            for k in range(len(new_list_r[i])):
                if math.sqrt((sample_g[i][j][0] - new_list_r[i][k][0])**2 + (sample_g[i][j][1] - new_list_r[i][k][1])**2) <= 0.0009:
                    if new_list_r[i][k] in R:
                        count += 1
                    G.append(sample_g[i][j])
                    R.append(new_list_r[i][k])
        print("Length of matched sample is ", len(G))
        print("Number of duplicate match sample in each pairs is", count)
        text3 = "The plot shows noise filtered sample in G & R " + str(number[i])
        plot_sample(G, R, text3)
    return all_delta_ra_dec


def data_shifting(r_data, delta):
     # --------------------------------all data shift section------------------------------------------
    for i in range(len(r_data)):
        for j in range(len(r_data[i])):
            r_data[i][j][0] = r_data[i][j][0] + delta[i][0]
            r_data[i][j][1] = r_data[i][j][1] + delta[i][1]


def data_rotate(g_data, r_data, center_point, theta):
    print("*****Rotating Section*****")
    number = [22,23,24,31,32,33]
    G, R = [], []
    for i in range(len(r_data)):
        # visualize before rotate
        text = "The plot shows data before rotate in G&R " + str(number[i])
        plot_sample(g_data[i], r_data[i], text)
        print('Rotating Sample R', number[i])
        # change all points to origin of picture center
        for j in range(len(r_data[i])):
            r_data[i][j][0] -= center_point[i][0]
            r_data[i][j][1] -= center_point[i][1]
            x = r_data[i][j][0]
            y = r_data[i][j][1]
        # convert the cartesian coordinate to polar coordinate (r,theta)
            r_data[i][j][0] = (x ** 2 + y ** 2) ** .5
            r_data[i][j][1] = math.degrees(math.atan2(y, x))
        # change parameter >> Theta <<
            r_data[i][j][1] += theta
            r = r_data[i][j][0]
            t = r_data[i][j][1]
        # convert polar coordinate to cartesian coordination
            r_data[i][j][0] = r*math.cos(math.radians(t))
            r_data[i][j][1] = r*math.sin(math.radians(t))
        # back to original
        for k in range(len(r_data[i])):
            r_data[i][k][0] += center_point[i][0]
            r_data[i][k][1] += center_point[i][1]
        # visualize after rotate the image
        text2 = "The plot shows rotated data in G & R " + str(number[i])
        plot_sample(g_data[i], r_data[i], text2)

        # data filtering for noiseless
        temp_G = []
        temp_R = []
        count = 0
        for m in range(len(g_data[i])):
            for n in range(len(r_data[i])):
                if math.sqrt((g_data[i][m][0] - r_data[i][n][0])**2 + (g_data[i][m][1] - r_data[i][n][1])**2) <= 0.00110:
                    if r_data[i][n] in temp_R:
                        count += 1
                    temp_G.append(g_data[i][m])
                    temp_R.append(r_data[i][n])
        G.append(temp_G)
        R.append(temp_R)
        print("Length of matched data rotation is ", len(G[i]))
        print("Number of duplicate match data in each pairs is", count)
        # print(G[i])
        text3 = "The plot shows filtered data in G & R " + str(number[i])
        plot_sample(G[i], R[i], text3)
    return G


def matching(key, data_G):
    # print("The number of key is: ",len(key))
    # print(key)
    # print(data_G)
    count = 0
    duplicate_check = []
    check = []
    matched_list = []
    for i in range(len(key)):
        for j in range(len(data_G)):
            for k in range(len(data_G[j])): 
                if(abs(key[i][0] - data_G[j][k][0]) <= 0.003 and abs(key[i][1] - data_G[j][k][1]) <= 0.003):
                    temp = []
                    temp1 = []
                    temp.append(float(data_G[j][k][0]))  # RA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                    temp.append(float(data_G[j][k][1]))  # DEC
                    temp.append(float(key[i][0])) #R.A. key
                    temp.append(float(key[i][1])) # Dec. key
                    temp.append(float(data_G[j][k][2]))  # MAG_APER
                    temp.append(data_G[j][k][3])  # filename
                    temp.append(data_G[j][k][4])  # object_number
                    matched_list.append(temp)
                    temp1.append(key[i][0])
                    temp1.append(key[i][1])
                    if(temp1 in duplicate_check):
                        count += 1
                        check.append(temp1)
                    else:
                        duplicate_check.append(temp1)
    print("Length of key is ", len(key))
    print("Length of matched list is ", len(matched_list))
    print("Length of duplicate matched list is ", count)
    print(check) 
    return matched_list


# splitter for split datasest in format of 'gX.txt' file only.
def splitter(g_file):
    splited = []
    myfile = open(g_file)
    info = myfile.readlines()
    for text in info:
        if(text[0] != '#'):
            text = text.rstrip('\n')
            text = re.findall(r'\S+', text)
            text.append(g_file)
            splited.append(text)
            myfile.close()
    return splited


def csv_key_export(key_list):
    with open('key_02.csv', 'w', newline='') as fp:
        a = csv.writer(fp, delimiter=',')
        data = [['RA', 'DEC', 'Vmag']]
        a.writerows(data)
        for i in range(len(key_list)):
            data = [[key_list[i][0], key_list[i][1], key_list[i][2]]]
            a.writerows(data)


def csv_export(matched_list):
    with open('matchData.csv', 'w', newline='') as fp:
        a = csv.writer(fp, delimiter=',')
        data = [['NUMBER', 'MAG_APER', 'MAGERR_APER', 'MAG_AUTO', 'MAGERR_AUTO', 'MAG_BEST', 'MAGERR_BEST', 'KRON_RADIUS', 'BACKGROUND', 'THRESHOLD', 'ISOAREA_IMAGE', 'ISOAREAF_IMAGE', 'ALPHAPEAK_J2000',
                 'DELTAPEAK_J2000', 'X_IMAGE', 'Y_IMAGE', 'FWHM_IMAGE', 'FWHM_WORLD', 'ELONGATION', 'ELLIPTICITY', 'CLASS_STAR', 'FLUX_RADIUS', 'FILE_NAME']]
        a.writerows(data)
        for i in range(len(matched_list)):
            m = matched_list[i]
            # use splitter function to split data in file name
            temp = splitter(m[5])
            data = [[temp[int(m[6])][0], temp[int(m[6])][1], temp[int(m[6])][2], temp[int(m[6])][3], temp[int(m[6])][4], temp[int(m[6])][5], temp[int(m[6])][6], temp[int(m[6])][7], temp[int(m[6])][8], temp[int(
                m[6])][9], temp[int(m[6])][10], temp[int(m[6])][11], temp[int(m[6])][12], temp[int(m[6])][13], temp[int(m[6])][14], temp[int(m[6])][15], temp[int(m[6])][16], temp[int(m[6])][17], temp[int(m[6])][18],
                temp[int(m[6])][19], temp[int(m[6])][20], temp[int(m[6])][21], m[5]]]
            a.writerows(data)


def histogram(matched_list):
    temp = []
    for i in range(len(matched_list)):
        temp.append(matched_list[i][2])
    plt.hist(temp, bins=20, ec='black')
    plt.title(
        'Histogram show visualization of luminosity function')
    plt.show()


def plot_sample(sample_g, sample_r, text):
    x_g, y_g, x_r, y_r = [], [], [], []
    for i in range(len(sample_g)):
        x_g.append(sample_g[i][0])
        y_g.append(sample_g[i][1])
    for j in range(len(sample_r)):
        x_r.append(sample_r[j][0])
        y_r.append(sample_r[j][1])
    plt.figure(figsize=(4.3, 10))

    plt.scatter(x_r, y_r, label='Sample R', c='white',
                edgecolors='crimson', marker='o', s=5)
    plt.scatter(x_g, y_g, label='Sample G', c='forestgreen', marker='.', s=5)
    plt.title(text)
    plt.xlabel('Alpha', fontsize=12)
    plt.ylabel('Delta', fontsize=12)
    plt.legend(fontsize=12)
    plt.show()


def plot_key_data(matched, key, text):
    x_matched, y_matched, x_key, y_key = [], [], [], []
    for i in range(len(matched)):
        for j in range(len(matched[i])):
            x_matched.append(matched[i][j][0])
            y_matched.append(matched[i][j][1])

    for j in range(len(key)):
        x_key.append(key[j][0])
        y_key.append(key[j][1])
    plt.figure(figsize=(4.3, 10))
    plt.scatter(x_key, y_key, label='Sample R', c='white',
                edgecolors='orangered', marker='o', s=5)
    plt.scatter(x_matched, y_matched, label='Sample G',
                c='forestgreen', marker='.', s=10)
    plt.title(text)
    plt.show()

def plot_filnal(matched, text):
    x_matched, y_matched, x_key, y_key = [], [], [], []
    for i in range(len(matched)):
        x_matched.append(matched[i][0])
        y_matched.append(matched[i][1])
        x_key.append(matched[i][2])
        y_key.append(matched[i][3])
    plt.figure(figsize=(4.3, 10))
    plt.scatter(x_key, y_key, label='Sample R', c='white',
                edgecolors='orangered', marker='o', s=5)
    plt.scatter(x_matched, y_matched, label='Sample G',
                c='forestgreen', marker='.', s=10)
    plt.title(text)
    plt.show()

if __name__ == "__main__":
    main()
