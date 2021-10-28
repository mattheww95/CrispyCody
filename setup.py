import setuptools
from setuptools import find_packages

setuptools.setup(
      name='CrispyCody',
      version=1.0,
      packages=find_packages(),
      description="Breadth of Coverage and depth calculato",
      scripts = ["src/CrispyCody"],
      zip_safe = False
      )
