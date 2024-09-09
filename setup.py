from setuptools import setup, find_packages

# 项目的元数据
NAME = "DongRadar"
VERSION = "0.0.1"
AUTHOR = "Dong Hongchang, Currently working at the Dingxi Meteorological Bureau"
AUTHOR_EMAIL = "dsg327@163.com"
DESCRIPTION = "weather radar tool"
URL = "https://github.com/dsg327/DongRadar.git"
LICENSE = "MIT"

# 项目依赖
REQUIRED = [
    "numpy",
    "struct",
    "pathlib",
]

# 读取 README 文件作为长描述
with open("README.md", "r", encoding="utf-8") as f:
    LONG_DESCRIPTION = f.read()

setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url=URL,
    license=LICENSE,
    packages=find_packages(),
    install_requires=REQUIRED,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
