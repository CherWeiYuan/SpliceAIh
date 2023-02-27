from setuptools import setup, find_packages

LONG_DESCRIPTION = \
"""
This programme runs SpliceAI on haplotypes
"""

setup(name="spliceaih",
      author="CHER_WEI_YUAN",
      version="1.0.0",
      author_email="E0031403@U.NUS.EDU",
      packages=find_packages(include=["spliceaih"]),
      entry_points={"console_scripts": ["spliceaih=spliceaih.__main__:main"]},
      package_data={"spliceaih": ["annotations/grch38.txt",
      "models/spliceai1.h5", "models/spliceai2.h5",
      "models/spliceai3.h5", "models/spliceai4.h5",
      "models/spliceai5.h5"]},
      url="https://github.com/CherWeiYuan/SpliceAIh/",
      description=("Bioinformatics commandline tool"),
      long_description=(LONG_DESCRIPTION))
