import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pathway-assessor-annapamma",
    version="0.0.1",
    author="Anna Calinawan",
    author_email="anna.calinawan@mssm.edu",
    description="For assessing the overrepresentation "
                "and underrepresentation of pathway genes "
                "in a given expression table",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/annapamma/PathwayAssessor",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
