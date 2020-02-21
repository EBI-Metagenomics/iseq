def get_version():
    import pkg_resources
    import re
    from os.path import realpath, dirname, join

    if __name__ == "__main__":
        filepath = join(dirname(realpath(__file__)), "..", "__init__.py")
        with open(filepath, "r", encoding="utf8") as f:
            content = f.read()
    else:
        content = pkg_resources.resource_string(__name__.split(".")[0], "__init__.py")
        content = content.decode()

    c = re.compile(r"__version__ *= *('[^']+'|\"[^\"]+\")")
    m = c.search(content)
    if m is None:
        return "unknown"
    return m.groups()[0][1:-1]
