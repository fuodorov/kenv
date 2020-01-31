from setuptools import setup, find_packages
import kenv

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='kenv',
    version=kenv.__version__,
    author='Vyacheslav Fedorov',
    author_email='fuodorov1998@gmail.com',
    description=kenv.__doc__,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/fuodorov/kenv',
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
