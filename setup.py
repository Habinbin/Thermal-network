from setuptools import setup, find_packages

setup(
    name="thermal_simulation",
    version="0.3",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A library for thermal simulation calculations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Habinbin/Thermal-network",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
