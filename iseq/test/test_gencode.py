import pytest

from iseq.gencode import GeneticCode


def test_genetic_code():
    with pytest.raises(ValueError):
        GeneticCode("Nao sei")

    with pytest.raises(ValueError):
        GeneticCode("Nao sei", "Sei la")

    gcode = GeneticCode("Invertebrate Mitochondrial")
    assert gcode.name == "Invertebrate Mitochondrial"
    assert gcode.alt_name == "SGC4"
    assert gcode.id == 5

    gcode = GeneticCode(alt_name="SGC4")
    assert gcode.name == "Invertebrate Mitochondrial"
    assert gcode.alt_name == "SGC4"
    assert gcode.id == 5

    gcode = GeneticCode(id=5)
    assert gcode.name == "Invertebrate Mitochondrial"
    assert gcode.alt_name == "SGC4"
    assert gcode.id == 5

    with pytest.raises(ValueError):
        GeneticCode(alt_name="")

    gcode = GeneticCode(name="Cephalodiscidae Mitochondrial")
    assert gcode.name == "Cephalodiscidae Mitochondrial"
    assert gcode.alt_name == ""
    assert gcode.id == 33
