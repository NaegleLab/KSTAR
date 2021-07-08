from setuptools import setup, find_packages


# Run setup
setup(
    name="kstar",
    version="0.1",
    author="Naegle Lab",
    author_email="kmn4mj@virginia.edu",
    url="https://github.com/NaegleLab/KSTAR",
    install_requires=['pandas', 'numpy', 'scipy', 'matplotlib', 'seaborn', 'statsmodels', 'biopython'],
    license='GNU General Public License v3',
    description="""KSTAR: Kinase-Substrate Transfer to Activity Relationships',
    long_description="KSTAR is an open-source software for estimating kinase activities from phosphoproteomic data. 
    KSTAR implements statistical and graph-theoretic approaches to produce a robust activity score that increases with 
    increasing evidence from a kinase's network.""",
    project_urls = {'Documentation': 'https://naeglelab.github.io/KSTAR/'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    python_requires=">=3.6",
    zip_safe = False
)

