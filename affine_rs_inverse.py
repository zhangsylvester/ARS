import matplotlib.pyplot as plt
import matplotlib.patches as patches
from tabloid import *
from affine_perm import *


def skew_to_list(T):
    """Return column index list for standard skew tableau T."""
    n = T.size()
    result = [0]*n
    for r, row in enumerate(T.tableau):
        for c, v in enumerate(row):
            if v is not None:
                result[v-1] = c+1
    return result

class Box:
    def __init__(self, top=None, bottom=None, left=None, right=None):
        self._top, self._bottom, self._left, self._right = top, bottom, left, right
    def top(self): return self._top
    def bottom(self): return self._bottom
    def left(self): return self._left
    def right(self): return self._right
    def set_top(self, v): self._top = v
    def set_bottom(self, v): self._bottom = v
    def set_left(self, v): self._left = v
    def set_right(self, v): self._right = v

class Window:
    def __init__(self, left_boundary, bottom_boundary):
        n = len(left_boundary)
        assert len(bottom_boundary)==n, "Boundaries length mismatch"
        self.n = n
        self._boxes = [[Box() for _ in range(n)] for _ in range(n)]
        self._dots = []
        for i in range(n):
            self._boxes[i][0].set_left(left_boundary[i])
            self._boxes[0][i].set_bottom(bottom_boundary[i])
        self._fill_edges()
        
    def _fill_edges(self):
        n = self.n
        for i in range(n):
            for j in range(n):
                b = self._boxes[i][j].bottom()
                l = self._boxes[i][j].left()
                if b!=l:
                    t, rt = b, l
                elif b==0:
                    t, rt = 0, 0
                else:
                    t, rt = b-1, b-1
                    if b==1:
                        self._dots.append((i,j))
                self._boxes[i][j].set_top(t)
                self._boxes[i][j].set_right(rt)
                if i+1<n: self._boxes[i+1][j].set_bottom(t)
                if j+1<n: self._boxes[i][j+1].set_left(rt)

    def boxes(self): return self._boxes
    def dots(self): return list(self._dots)
    def top_boundary(self): return [self._boxes[self.n-1][j].top() for j in range(self.n)]
    def right_boundary(self): return [self._boxes[i][self.n-1].right() for i in range(self.n)]
    def bottom_boundary(self): return [self._boxes[0][j].bottom() for j in range(self.n)]
    def left_boundary(self): return [self._boxes[i][0].left() for i in range(self.n)]

    def display(self):
        """
        Display this window's edge labels and dots with a dashed grid.
        """
        fig, ax = plt.subplots()
        n = self.n

        # Draw dashed grid lines
        for r in range(n+1):
            ax.plot([0, n], [r, r], linestyle='--', color='gray')
        for c in range(n+1):
            ax.plot([c, c], [0, n], linestyle='--', color='gray')

        eps = 0.2
        # Annotate each box's edges
        for i in range(n):
            for j in range(n):
                box = self._boxes[i][j]
                # Bottom edge
                b = box.bottom()
                if b is not None:
                    ax.text(j + 0.5, i + eps, str(b), ha='center', va='bottom')
                # Top edge
                t = box.top()
                if t is not None:
                    ax.text(j + 0.5, i + 1 - eps, str(t), ha='center', va='top')
                # Left edge
                l = box.left()
                if l is not None:
                    ax.text(j + eps, i + 0.5, str(l), ha='left', va='center')
                # Right edge
                rt = box.right()
                if rt is not None:
                    ax.text(j + 1 - eps, i + 0.5, str(rt), ha='right', va='center')

        # Plot dots
        for (i, j) in self._dots:
            ax.plot(j + 0.5, i + 0.5, 'ko')

        ax.set_aspect('equal')
        ax.axis('off')
        plt.show()


class InverseGrowthDiagram:
    def __init__(self, P, Q, lam, k):
        # 1) P, Q must be Tabloid
        if not isinstance(P, Tabloid) or not isinstance(Q, Tabloid):
            raise TypeError("P and Q must be Tabloids.")
        mu = P.shape()
        # 2) matching shapes
        if Q.shape() != mu:
            raise ValueError("P and Q must have the same shape.")
        # 3) skew tableaux standard
        skewP = standard_skew_tableau(lam, P.rows)
        skewQ = standard_skew_tableau(lam, Q.rows)
        if not skewP.is_standard():
            raise ValueError("lambda+P must be standard.")
        if not skewQ.is_standard():
            raise ValueError("lambda+Q must be standard.")
        # 4) length condition
        if len(lam) > len(mu):
            raise ValueError("lambda cannot have more parts than mu.")
        # 5) strict exceed or extra rows
        if not (len(mu) > len(lam) or any(mu[k] > lam[k] for k in range(min(len(mu), len(lam))))):
            raise ValueError("lambda <= mu.")
        self.n = sum(mu)
        self.P, self.Q, self.lam, self.k, = P, Q, lam, k
        # Prepare P_list and Q_list
        self.P_list = skew_to_list(standard_skew_tableau(lam, P.rows))
        self.Q_list = skew_to_list(standard_skew_tableau(lam, Q.rows))
        # Build successive windows
        self._windows = []
        current = Window(self.Q_list[::-1], self.P_list[::-1])
        while True:
            self._windows.append(current)
            if all(v == 0 for v in current.top_boundary()) and all(v == 0 for v in current.right_boundary()):
                break
            current = Window(current.right_boundary(), current.top_boundary())
        # Aggregate edge boxes into a 2D list: n rows, m*n columns
        n = self._windows[0].n
        m = len(self._windows)
        self._edge_boxes = []
        for i in range(n):
            row_boxes = []
            for w in self._windows:
                row_boxes.extend(w.boxes()[i])
            self._edge_boxes.append(row_boxes)
        # Aggregate dot positions with window offset
        self._dots = []
        for k, w in enumerate(self._windows):
            for (r, c) in w.dots():
                self._dots.append((r, c + k * n))

    def windows(self):
        return list(self._windows)

    def edge_labels(self):
        """Return a 2D list of Box objects: n rows by m*n columns."""
        return [list(row) for row in self._edge_boxes]

    def dots(self):
        """Return list of (row, col) positions of dots in the aggregated grid."""
        return list(self._dots)

    # def display_edge_labels(self):
    #     """Visualize the full growth diagram edge labels and dots,
    #     with dashed grid, edges labeled at midpoints, omitting repeats,
    #     and rotated counterclockwise 90°."""
    #     # Original grid: n rows x total_cols columns
    #     n = len(self._edge_boxes)
    #     total_cols = len(self._edge_boxes[0]) if n>0 else 0
    #     # Build rotated grid: rows'=total_cols, cols'=n
    #     rotated = [ [self._edge_boxes[j][total_cols-1-i] for j in range(n)]
    #                 for i in range(total_cols) ]
    #     rows_rot, cols_rot = total_cols, n
    #     fig, ax = plt.subplots()
    #     # Draw dashed grid
    #     for r in range(rows_rot+1):
    #         ax.plot([0, cols_rot], [r, r], linestyle='--', color='gray')
    #     for c in range(cols_rot+1):
    #         ax.plot([c, c], [0, rows_rot], linestyle='--', color='gray')
    #     # Track drawn edge labels to omit repeats
    #     drawn = set()
    #     eps = 0.2
    #     # Draw edge labels
    #     for r in range(rows_rot):
    #         for c in range(cols_rot):
    #             box = rotated[r][c]
    #             # bottom edge of rotated corresponds to original left
    #             label = box.bottom()
    #             seg = (r, c, 'h')
    #             if label is not None and (seg, label) not in drawn:
    #                 ax.text(c+0.5, r+eps, str(label), ha='center', va='bottom')
    #                 drawn.add((seg, label))
    #             # top edge -> original right
    #             label = box.top()
    #             seg = (r+1, c, 'h')
    #             if label is not None and (seg, label) not in drawn:
    #                 ax.text(c+0.5, r+1-eps, str(label), ha='center', va='top')
    #                 drawn.add((seg, label))
    #             # left edge -> original bottom
    #             label = box.left()
    #             seg = (r, c, 'v')
    #             if label is not None and (seg, label) not in drawn:
    #                 ax.text(c+eps, r+0.5, str(label), ha='left', va='center')
    #                 drawn.add((seg, label))
    #             # right edge -> original top
    #             label = box.right()
    #             seg = (r, c+1, 'v')
    #             if label is not None and (seg, label) not in drawn:
    #                 ax.text(c+1-eps, r+0.5, str(label), ha='right', va='center')
    #                 drawn.add((seg, label))
    #     # Plot dots (rotate positions)
    #     for (r0, c0) in self._dots:
    #         # original at (r0,c0), window offset applied
    #         # in rotated: new_r = total_cols-1-c0, new_c = r0
    #         new_r = total_cols-1 - c0
    #         new_c = r0
    #         ax.plot(new_c+0.5, new_r+0.5, 'ko')
    #     ax.set_aspect('equal')
    #     ax.axis('off')
    #     plt.show()
    def inverse_rs(self):
        h = len(self.windows())*self.n
        perm = []
        for i in range(self.n):
            for dot in self.dots():
                if dot[0] == i:
                    perm += [h - dot[1]]
        fake_aff_perm = affine_permutation(perm)
        aff_index = (self.k - 1)*self.n - sum(self.lam)
        offset = aff_index - fake_aff_perm.index()
        return fake_aff_perm.shift(offset)
    def display_edge_labels(self):
        """
        Display the full growth diagram edge labels and dots
        without rotation or duplicate suppression.
        """
        fig, ax = plt.subplots()
        n = len(self._edge_boxes)
        total_cols = len(self._edge_boxes[0]) if n > 0 else 0
        # Draw dashed grid
        for r in range(n+1):
            ax.plot([0, total_cols], [r, r], linestyle='--', color='gray')
        for c in range(total_cols+1):
            ax.plot([c, c], [0, n], linestyle='--', color='gray')
        eps = 0.2
        # Annotate edges for every box
        for i in range(n):
            for j in range(total_cols):
                box = self._edge_boxes[i][j]
                # bottom edge
                val = box.bottom()
                if val is not None:
                    ax.text(j+0.5, i+eps, str(val), ha='center', va='bottom')
                # top edge
                val = box.top()
                if val is not None:
                    ax.text(j+0.5, i+1-eps, str(val), ha='center', va='top')
                # left edge
                val = box.left()
                if val is not None:
                    ax.text(j+eps, i+0.5, str(val), ha='left', va='center')
                # right edge
                val = box.right()
                if val is not None:
                    ax.text(j+1-eps, i+0.5, str(val), ha='right', va='center')
        # Plot dots
        for (i, j) in self._dots:
            ax.plot(j+0.5, i+0.5, 'ko')
        ax.set_aspect('equal')
        ax.axis('off')
        plt.show()