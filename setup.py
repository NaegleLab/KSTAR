from setuptools import setup, find_packages


# Run setup
setup(
    name="kstar",
    version="1.0.0",
    author="Naegle Lab",
    author_email="kmn4mj@virginia.edu",
    url="https://github.com/NaegleLab/KSTAR",
    install_requires=['pandas==2.1.*', 'numpy==1.26.*', 'scipy==1.11.*', 'matplotlib==3.8.*', 'seaborn==0.13.*', 'biopython==1.81.*','requests==2.31.*'],
    license='GNU General Public License v3',
    description='KSTAR: Kinase-Substrate Transfer to Activity Relationships',
    long_description="""KSTAR is an open-source software for estimating kinase activities from phosphoproteomic data. 
    KSTAR implements statistical and graph-theoretic approaches to produce a robust activity score that increases with increasing evidence from a kinase's network.""",
    project_urls = {'Documentation': 'https://naeglelab.github.io/KSTAR/'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    include_package_data = True,
    python_requires=">=3.11",
    zip_safe = False
)
