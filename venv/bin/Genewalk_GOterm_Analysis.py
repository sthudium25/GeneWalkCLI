from GenewalkObj import GenewalkObj
from typing import Dict


def main():
    print("Enter project name: ")
    project_title: str = input()
    gw = GenewalkObj(project_name=project_title)
    gw.load_results()


if __name__ == "__main__":
    main()
