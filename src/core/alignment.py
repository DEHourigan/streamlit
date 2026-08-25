"""BLAST alignment display helpers."""

def format_alignment(qseq: str, sseq: str, qstart: int, sstart: int, width: int = 60) -> str:
    blocks, q_position, s_position = [], int(qstart), int(sstart)
    for index in range(0, len(qseq), width):
        qblock, sblock = qseq[index:index + width], sseq[index:index + width]
        matches = "".join(" " if "-" in (q, s) else "│" if q == s else "·" for q, s in zip(qblock, sblock))
        q_end = q_position + len(qblock.replace("-", "")) - 1
        s_end = s_position + len(sblock.replace("-", "")) - 1
        blocks.append(f"Query  {q_position:<6} {qblock}  {q_end}\n               {matches}\nHit    {s_position:<6} {sblock}  {s_end}")
        q_position, s_position = q_end + 1, s_end + 1
    return "\n\n".join(blocks)
