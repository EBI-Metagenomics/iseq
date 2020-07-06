import pytest

from iseq.gencode import GeneticCode, get_genetic_code


def test_genetic_code():
    with pytest.raises(KeyError):
        get_genetic_code("Nao sei")

    with pytest.raises(ValueError):
        get_genetic_code("Nao sei", "Sei la")

    gcode = get_genetic_code("Invertebrate Mitochondrial")
    assert gcode.name == "Invertebrate Mitochondrial"
    assert gcode.alt_name == "SGC4"
    assert gcode.id == 5

    gcode = get_genetic_code(alt_name="SGC4")
    assert gcode.name == "Invertebrate Mitochondrial"
    assert gcode.alt_name == "SGC4"
    assert gcode.id == 5

    gcode = get_genetic_code(id=5)
    assert gcode.name == "Invertebrate Mitochondrial"
    assert gcode.alt_name == "SGC4"
    assert gcode.id == 5

    with pytest.raises(KeyError):
        get_genetic_code(alt_name="")

    gcode = get_genetic_code(name="Cephalodiscidae Mitochondrial")
    assert gcode.name == "Cephalodiscidae Mitochondrial"
    assert gcode.alt_name == ""
    assert gcode.id == 33
