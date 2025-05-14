class affine_permutation:
    def __init__(self, entries):
        """
        entries: list of integers of length n
        Requires that { entry mod n } = {0,1,...,n-1}
        """
        self.entries = list(entries)
        self.n = len(self.entries)
        mods = sorted([x % self.n for x in self.entries])
        if mods != list(range(self.n)):
            raise ValueError(
                f"Entries mod {self.n} must be a permutation of 0..{self.n-1}, got {mods}"
            )

    def __repr__(self):
        return f"affine_permutation {self.entries} "

    def __getitem__(self, i):
        """Return the i-th entry (0-based)."""
        return self.entries[i]

    def __len__(self):
        return self.n

    def index(self):
        """
        Compute the affine index:
          i = (sum(entries) - n*(n+1)/2) / n
        Returns an integer.
        """
        total = sum(self.entries)
        offset = self.n * (self.n + 1) // 2
        return (total - offset) // self.n

    def shift(self, x):
        """
        Return a new affine_permutation with each entry increased by x.
        """
        if not isinstance(x, int):
            raise TypeError("Shift amount must be an integer.")
        return affine_permutation([e + x for e in self.entries])
    
    def rotate(self, m):
        """
        Rotate the affine permutation:
          if m == 0: return self
          if m > 0: apply forward rotation m times: [w1,...,w_n] -> [w2-1,...,w_n-1,w1+n-1]
          if m < 0: apply backward rotation |m| times: [w1,...,w_n] -> [w_n-n+1,w1+1,...,w_{n-1}+1]
        Returns a new affine_permutation.
        """
        if not isinstance(m, int):
            raise TypeError("Rotation count must be an integer.")
        # no rotation
        if m == 0:
            return affine_permutation(self.entries)
        res = list(self.entries)
        # forward rotations
        if m > 0:
            for _ in range(m % self.n):
                old = res[:]
                res = [old[i+1] - 1 for i in range(self.n - 1)]
                res.append(old[0] + self.n - 1)
        else:
            # backward rotations
            for _ in range((-m) % self.n):
                old = res[:]
                # new first = w_n - (n-1)
                first = old[-1] - (self.n - 1)
                rest = [old[i] + 1 for i in range(self.n - 1)]
                res = [first] + rest
        return affine_permutation(res)