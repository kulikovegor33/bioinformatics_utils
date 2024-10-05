from typing import List, Union, Dict, Tuple

from scripts.script_dna_rna_tools import is_valid_sequence, transcribe, reverse
from scripts.script_dna_rna_tools import complement, reverse_complement
from scripts.script_fastq_filter import gc_content, average_quality


def run_dna_rna_tools(
    *args: Union[str, List[str]]
) -> Union[str, List[str], None]:

    """
    Выполняет указанные процедуры (транскрипция,
реверс, комплементарность) для заданных
последовательностей ДНК/РНК.

    Args:
        *args: Последовательности ДНК/РНК и
процедура, которую нужно выполнить.
               Последний аргумент должен быть строкой,
обозначающей процедуру.

    Returns:
        Результат выполнения процедуры:
        - Если процедура "transcribe", возвращает
транскрибированную последовательность.
        - Если процедура "reverse",
возвращает развёрнутую последовательность.
        - Если процедура "complement",
возвращает комплементарную последовательность.
        - Если процедура "reverse_complement",
возвращает обратную комплементарную последовательность.

        Если последовательностей одна, возвращает
строку, если больше одной - возвращает список строк.
        Если последовательностей нет, возвращает
None.

    Raises:
        ValueError: Если последовательность
недействительна или процедура неизвестна.
    """

    if not args:
        return None

    procedure = args[-1]
    sequences = args[:-1]

    for seq in sequences:
        if not is_valid_sequence(seq):
            raise ValueError(f"Invalid sequence: {seq}")

    if procedure == "transcribe":
        result = [transcribe(seq) for seq in sequences]
    elif procedure == "reverse":
        result = [reverse(seq) for seq in sequences]
    elif procedure == "complement":
        result = [complement(seq) for seq in sequences]
    elif procedure == "reverse_complement":
        result = [reverse_complement(seq) for seq in sequences]
    else:
        raise ValueError(f"Unknown procedure: {procedure}")

    return result if len(result) > 1 else result[0]


def filter_fastq(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Tuple[float, float] = (0.0, 100.0),
    length_bounds: Tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0
) -> Dict[str, Tuple[str, str]]:
    """
    Фильтрует последовательности в формате FASTQ
на основе заданных критериев.

    Args:
        seqs: Словарь с именами последовательностей
в качестве ключей и кортежами (последовательность,
качество) в качестве значений.
              Структура: {имя_последовательности:
(последовательность, качество)}
        gc_bounds: Интервал GC-состава (в процентах)
для фильтрации.
                   По умолчанию (0, 100) — все риды
сохраняются.
                   Если передать одно число,
считается верхняя граница.
                   Примеры: (20, 80) — сохраняем
риды с GC составом от 20 до 80%,
                   (44.4) — сохраняем риды с GC
составом меньше 44.4%.
        length_bounds: Интервал длины для
фильтрации. По умолчанию (0, 2**32).
                       Аналогично gc_bounds.
        quality_threshold: Пороговое значение
среднего качества рида для фильтрации.
                           По умолчанию равно 0
(шкала phred33). Риды со средним качеством ниже
порога отбрасываются.

    Returns:
        Словарь с отфильтрованными
последовательностями в формате {имя_последовательности: (последовательность,
качество)}.
    """

    filtered_seqs = {}

    # Обработка gc_bounds
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0.0, float(gc_bounds))

    # Обработка length_bounds
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    for name, (sequence, quality) in seqs.items():
        gc_perc = gc_content(sequence)
        seq_length = len(sequence)
        avg_quality = average_quality(quality)

        if (gc_bounds[0] <= gc_perc <= gc_bounds[1] and
                length_bounds[0] <= seq_length <= length_bounds[1] and
                avg_quality >= quality_threshold):
            filtered_seqs[name] = (sequence, quality)

    return filtered_seqs
