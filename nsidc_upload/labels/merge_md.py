import csv
import os
import glob
import argparse
import numpy as np


def merge_csv_files(src_dir, dst_dir):
    '''
    Takes a list of csv files (full file path) and merges them into a
    single output.
    '''
    # Create a list of all the md files in the given directory
    csv_list = glob.glob(os.path.join(src_dir, '*_md.csv'))
    csv_list.sort()     # This will sort by frame number


    # extract the flight name (i.e. date) from one of the csv files
    fname = os.path.basename(csv_list[0])
    flight_id = fname.rsplit('_', 2)[0].split('/')[-1]
    yyyy, mm, dd = flight_id.split('_')

    aggregate_data = []

    # Read each file and add it to the aggregate list
    for file in csv_list:
        data = read_md_file(file)
        aggregate_data.append(data)

    dst_name = 'GR_{}{}{}_metadata.csv'.format(yyyy, mm, dd)
    dst_file = os.path.join(dst_dir, dst_name)

    # Write the aggregated data to a single output file
    write_merged_md(dst_file, aggregate_data)


def read_md_file(md_file):
    '''
    Reads a single metadata file and returns the relevant information
    '''
    qa = 0
    snow, gray, pond, ocean, shadow = 0, 0, 0, 0, 0

    # Get the frame number from the filename
    fname = os.path.basename(md_file)
    frame = int(fname.split('_')[3])

    with open(md_file, 'r') as md:
        csv_reader = csv.reader(md)
        # Remove the header
        try:
            next(csv_reader)
            for row in csv_reader:
                qa = float(row[0])
                snow = float(row[1])
                gray = float(row[2])
                pond = float(row[3])
                ocean = float(row[4])
                try:
                    shadow = float(row[5])
                except IndexError:
                    shadow = 0
        except Exception as e: # Make sure empty files don't crash the tool
            print('Caught exception: ' + str(e))
            print('File: ' + md_file)

    data = [frame, qa, snow, gray, pond, ocean, shadow]

    return data


def write_merged_md(dst_file, aggregate_data):

    print("Writing to {}...".format(dst_file))

    with open(dst_file, 'w') as md:
        csv_writer = csv.writer(md)
        csv_writer.writerow(["Frame", "Quality Score", "White Ice",
                             "Gray Ice", "Melt Ponds", "Open Water", "Shadow"])
        csv_writer.writerows(aggregate_data)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("src_dir",
                        help="Source directory containing *_md.csv files.")
    parser.add_argument("--dst_dir", default=None,
                        help="folder to place the merged data file")

    args = parser.parse_args()

    #src_dir = os.path.dirname(args.src_dir)
    src_dir = args.src_dir
    dst_dir = args.dst_dir

    if dst_dir is None:
        dst_dir = os.path.split(src_dir)[0]
    #else:
    #    dst_dir = os.path.dirname(dst_dir)

    merge_csv_files(src_dir, dst_dir)


if __name__ == '__main__':
    main()
