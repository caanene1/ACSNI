import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as rq:
    install_requires = rq.read()


setuptools.setup(
    name="ACSNI",
    version="1.0.0",
    scripts=["ACSNI", "ACSNI-derive", "ACSNI-get", "ACSNI-split"],
    author="Chinedu A. Anene",
    author_email="caanenedr@outlook.com",
    description="automatic context-specific network inference",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/caanene1/ACSNI",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    install_requires=install_requires,
)

# Build this package with >> python3 setup.py sdist bdist_wheel