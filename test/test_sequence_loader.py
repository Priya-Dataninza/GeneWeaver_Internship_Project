from src.sequence_loader import load_sequence

def test_valid_sequence():
    result = load_sequence("ATGC")
    assert result == "ATGC"
    print("Valid sequence test passed!")

def test_lowercase_sequence():
    result = load_sequence("atgc")
    assert result == "ATGC"
    print("Lowercase test passed!")

def test_empty_sequence():
    try:
        load_sequence("")
    except ValueError:
        print("Empty sequence test passed!")

def test_invalid_sequence():
    try:
        load_sequence("ATGX")
    except ValueError:
        print("Invalid sequence test passed!")

if __name__ == "__main__":
    test_valid_sequence()
    test_lowercase_sequence()
    test_empty_sequence()
    test_invalid_sequence()

    print("All tests completed!")

