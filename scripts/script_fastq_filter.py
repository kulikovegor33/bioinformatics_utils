def gc_content(sequence: str) -> float:
    """
    Вычисляет содержание гуанина и цитозина (GC) в последовательности.

    Args:
        sequence (str): Последовательность ДНК.

    Returns:
        float: Процентное содержание GC в последовательности.
    """
    g_count = sequence.count('G') + sequence.count('g')
    c_count = sequence.count('C') + sequence.count('c')
    total_count = len(sequence)

    if total_count == 0:
        return 0.0

    return (g_count + c_count) / total_count * 100


def average_quality(quality: str) -> float:
    """
    Вычисляет среднее качество последовательности на основе кодов ASCII.

    Args:
        quality (str): Строка качества (обычно в формате Phred).

    Returns:
        float: Среднее значение качества.
    """
    return sum(ord(q) - 33 for q in quality) / len(quality)
