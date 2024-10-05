VALID_DNA = set("ATGCatgc")
VALID_RNA = set("AUGCaugc")
COMPLEMENT_MAP_DNA = str.maketrans("ATGCatgc", "TACGtacg")
COMPLEMENT_MAP_RNA = str.maketrans("AUGCaugc", "UACGuacg")


def is_valid_sequence(sequence: str) -> bool:
    """
    Проверяет, является ли последовательность допустимой ДНК или РНК.

    Args:
        sequence (str): Последовательность для проверки.

    Returns:
        bool: True, если последовательность допустима, иначе False.
    """
    return all(base in VALID_DNA or base in VALID_RNA for base in sequence)


def transcribe(dna_sequence: str) -> str:
    """
    Транскрибирует ДНК в РНК, заменяя 'T' на 'U'.

    Args:
        dna_sequence (str): Последовательность ДНК.

    Returns:
        str: Транскрибированная последовательность РНК.
    """
    return dna_sequence.replace("T", "U").replace("t", "u")


def reverse(sequence: str) -> str:
    """
    Возвращает обратную последовательность.

    Args:
        sequence (str): Исходная последовательность.

    Returns:
        str: Обратная последовательность.
    """
    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Находит комплементарную последовательность.

    Args:
        sequence (str): Исходная последовательность.

    Returns:
        str: Комплементарная последовательность.
    """
    if "U" in sequence or "u" in sequence:
        return sequence.translate(COMPLEMENT_MAP_RNA)
    else:
        return sequence.translate(COMPLEMENT_MAP_DNA)


def reverse_complement(sequence: str) -> str:
    """
    Находит обратную комплементарную последовательность.

    Args:
        sequence (str): Исходная последовательность.

    Returns:
        str: Обратная комплементарная последовательность.
    """
    return reverse(complement(sequence))
