
class Tabloid:
    def __init__(self, rows):
        # rows: list of lists of integers
        self.rows = rows

    def display(self):
        for row in self.rows:
            print(' '.join(str(x) for x in row))

    def __repr__(self):
        return f"Tabloid({self.rows})"

    def shape(self):
        """Return the shape as list of row lengths."""
        return [len(row) for row in self.rows]

class skew_tableau:
    def __init__(self, mu, rows):
        """
        mu: list of nonnegative ints; rows: list of lists of values
        Builds each row with mu[i] leading None then rows[i].
        """
        self.mu = mu
        self.rows_input = rows
        self.tableau = []
        for i, r in enumerate(rows):
            if i < len(mu):
                self.tableau.append([None]*mu[i] + list(r))
            else:
                self.tableau.append(list(r))

    def content(self):
        """Flatten non-None entries."""
        return [e for row in self.tableau for e in row if e is not None]

    def rows(self):
        """Return rows without None."""
        return [[e for e in row if e is not None] for row in self.tableau]

    def columns(self):
        """Return columns without None."""
        maxc = max(len(row) for row in self.tableau)
        cols = []
        for j in range(maxc):
            col = []
            for row in self.tableau:
                if j < len(row) and row[j] is not None:
                    col.append(row[j])
            cols.append(col)
        return cols

    def is_standard(self):
        """Check consecutive content 1..n and rows/cols strictly increasing."""
        cont = sorted(self.content())
        n = len(cont)
        if cont != list(range(1, n+1)):
            return False
        for row in self.rows():
            if any(row[i] >= row[i+1] for i in range(len(row)-1)):
                return False
        for col in self.columns():
            if any(col[i] >= col[i+1] for i in range(len(col)-1)):
                return False
        return True

    def __str__(self):
        """Pretty-print with '*' for None."""
        maxc = max(len(row) for row in self.tableau)
        widths = [0]*maxc
        for j in range(maxc):
            entries = [str(row[j]) for row in self.tableau if j < len(row) and row[j] is not None]
            widths[j] = max((len(s) for s in entries), default=1)
        lines = []
        for row in self.tableau:
            parts = []
            for j in range(len(row)):
                if row[j] is not None:
                    parts.append(str(row[j]).rjust(widths[j]))
                else:
                    parts.append('*'*widths[j])
            lines.append(' '.join(parts))
        return '\n'.join(lines)


def tabloid_to_skew(T, mu):
    """Convert Tabloid to skew_tableau using its rows."""
    return skew_tableau(mu, T.rows)

class standard_skew_tableau(skew_tableau):
    def __init__(self, mu, rows):
        # enforce mu weakly decreasing
        if any(mu[k] < mu[k+1] for k in range(len(mu)-1)):
            raise ValueError("mu must be weakly decreasing.")
        super().__init__(mu, rows)
        if not self.is_standard():
            raise ValueError("Not a standard skew tableau.")

    def size(self):
        """Number of entries."""
        return len(self.content())

    def toggle(self, i):
        cont = self.content()
        if i not in cont or (i+1) not in cont:
            raise IndexError(f"{i} or {i+1} missing.")
        pos = {}
        for r, row in enumerate(self.tableau):
            for c, v in enumerate(row):
                if v == i or v == i+1:
                    pos[v] = (r, c)
        new_table = [list(r) for r in self.tableau]
        (r1, c1), (r2, c2) = pos[i], pos[i+1]
        new_table[r1][c1], new_table[r2][c2] = new_table[r2][c2], new_table[r1][c1]
        new_rows = [[e for e in row if e is not None] for row in new_table]
        try:
            return standard_skew_tableau(self.mu, new_rows)
        except ValueError:
            return self

    def Toggle(self, indices):
        result = self
        for k in indices:
            result = result.toggle(k)
        return result

    def q(self, i):
        n = self.size()
        if not (1 <= i < n):
            raise IndexError(f"i must be 1..{n-1}.")
        seq = list(range(1, i+1)) + list(range(i-1, 0, -1))
        return self.Toggle(seq)

    def evac(self, i, j):
        n = self.size()
        if not (1 <= i < j < n):
            raise IndexError("Require 1 <= i < j < size().")
        return self.q(j-1).q(j-i).q(j-1)