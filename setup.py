import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as rq:
    install_requires = rq.read()


setuptools.setup(
    name="ACSNI",
    version="1.0.4",
    scripts=["bin/ACSNI-run", "bin/ACSNI-derive", "bin/ACSNI-get", "bin/ACSNI-split"],
    author="Chinedu A. Anene",
    collaborator="Faraz Khan",
    author_email="caanenedr@outlook.com",
    description="automatic context-specific network inference",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/caanene1/ACSNI",
    download_url = "https://github.com/caanene1/ACSNI/releases/download/1.0.4/ACSNI-1.0.4-py3-none-any.whl",
    packages=setuptools.find_packages(include=["ACSNI"]),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    install_requires=install_requires,
    python_requires='>=3.6',
)

# Build >> python3 setup.py sdist bdist_wheel
# Upload >> sudo twine upload dist/*
