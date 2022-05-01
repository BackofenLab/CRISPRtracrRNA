import os
from os import path


def folder_maker(folder_output_path, folder_temp_path):
    """
    This function creates the output folder and the temporary folder.
    """

    # Check if the output folder exists. If not, create it.
    if not path.exists(folder_output_path):
        os.makedirs(folder_output_path)

    # Check if the temporary folder exists. If not, create it.
    if folder_temp_path:
        if not path.exists(folder_temp_path):
            os.makedirs(folder_temp_path)
