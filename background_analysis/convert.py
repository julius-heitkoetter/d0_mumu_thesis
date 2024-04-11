import os
import sys

target_path = "./"
maps = {"small-inclusive.root":
    [
        "InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM",
    ]
}


def main():
    """
    __main__
    """
    base_path = '/data/submit/juliush/d0measurement/FlatNtuples/529/dstar/'
    for target_file in maps:
        cmd = "hadd "+target_file+" "
        for input_files in maps[target_file]:
            cmd += base_path+input_files+"/*.root "
        os.system(cmd)

if __name__ == "__main__":
    sys.exit(main())