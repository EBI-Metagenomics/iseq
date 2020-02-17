from setuptools import setup

if __name__ == "__main__":
    console_scripts = ["iseq = iseq:cli"]
    setup(entry_points={"console_scripts": console_scripts},)
