import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pathway-assessor",
    version="0.0.2",
    author="Anna Calinawan",
    author_email="anna.calinawan@mssm.edu",
    description="For assessing the overrepresentation "
                "and underrepresentation of pathway genes "
                "in a given expression table",
    licence="MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/annapamma/PathwayAssessor",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
          "numpy",
          "scipy",
          "pandas"
    ],
    include_package_data=True,
)
