# Implement Hungarian algorithm

def solve_simple_assignment(quals):
  def compute_initial():
    # Find an initial assignment by simply using the first unassigned job in each row
    assignment = [[False for elem in row] for row in quals]
    assigned_cols = [False] * len(quals[0])
    for a_row, q_row in zip(assignment,quals):
      for i, q_elem in enumerate(q_row):
        if q_elem == True and not assigned_cols[i]:
          a_row[i] = True
          assigned_cols[i] = True
          break
    return assignment

  def enumerate_transfers(assignment):
    # Determine all transfers using BFS (NOT the DFS proposed in the Kuhn 1955)

    improving_transfers = []
    nonimproving_transfers = []
    incomplete_transfers = []
    
    # Find transfer starting locations
    for j in xrange(len(quals[0])):
      if all([not assignment[i][j] for i in xrange(len(quals))]):
        qual_rows = [i for i in range(len(quals)) if quals[i][j]]
        incomplete_transfers.extend([[(i,j)] for i in qual_rows])

    # Continue incomplete transfers until none are left
    while len(incomplete_transfers) > 0:
      new_incomplete_transfers = []
      for transfer in incomplete_transfers:
        # Search the last row in the transfer for an assigned column
        # If none, this row can be assigned to this column, so the transfer is complete (and an improvement)
        # If an assigned column is found, extend the transfer accordingly
        i,j = transfer[-1]
        if all([not a for a in assignment[i]]):
          improving_transfers.append(transfer)
          continue
        new_j = [new_j for new_j,a in enumerate(assignment[i]) if a][0]
        transfer.append((i,new_j))

        # Search this new column (new_j) for an unassigned row.
        # If none, the transfer is complete (but NOT an improvement)
        # For each unassigned row found, a new transfer to be extended is added to the growing list
        used_rows = set([i for i,_ in transfer])
        new_rows = [new_i for new_i in range(len(quals)) if quals[new_i][new_j] and new_i not in used_rows]
        if len(new_rows) == 0:
          nonimproving_transfers.append(transfer)
        else:
          new_incomplete_transfers.extend([transfer+[(new_i,new_j)] for new_i in new_rows])
      incomplete_transfers = new_incomplete_transfers

    return (improving_transfers, nonimproving_transfers)
          

  # Determine initial assignments
  assignment = compute_initial()

  # Keep searching for improving transfers (i.e. transfers w/ odd # entries),
  # improving the assignment each time one is found.
  improving_transfers, nonimproving_transfers = enumerate_transfers(assignment)
  while len(improving_transfers) > 0:
    transfer = improving_transfers[0]
    for i,j in transfer:
      assignment[i][j] = not assignment[i][j]

    improving_transfers, nonimproving_transfers = enumerate_transfers(assignment)

  # When none are found, mark essential rows (those found in at least one transfer)
  # and essential columns (columns assigned to nonessential rows)
  # and return a tuple: (the optimal assignment, essential rows, essential columns)
  essential_rows = [False] * len(quals)
  for t in nonimproving_transfers:
    for i,j in t:  essential_rows[i] = True

  essential_cols = [False] * len(quals[0])
  for i, row in enumerate(assignment):
    if essential_rows[i]:  continue
    for j, v in enumerate(row):
      if v:  essential_cols[j] = True
  
  return (assignment, essential_rows, essential_cols)

def solve_general_assignment(ratings):
  def compute_initial():
    a = [max(row) for row in ratings]
    b = [max([ratings[i][j] for i in range(len(ratings))]) for j in range(len(ratings[0]))]

    if sum(a) <= sum(b):
      return a, [0]*len(b)
    else:
      return [0]*len(a), b

  def solve_corresponding_simple_assignment(u, v):
    # Compute the qualification matrix for the corresponding simple assignment problem,
    # given the budgets assigned to rows (u) and columns (v)
    quals = []
    for i, row in enumerate(ratings):
      q_row = []
      for j, r_ij in enumerate(row):
        q_row.append(u[i]+v[j]==r_ij)
      quals.append(q_row)

    # Solve simple assignment problem corresponding to this qualification matrix
    return solve_simple_assignment(quals)

  # Compute initial budget
  u,v = compute_initial()

  def check_optimal_solution(assignment, essential_rows, essential_cols):
    return all(essential_rows) or all(essential_cols) ## not sure about this

  # Solve simple assignment problem until this solution solves the general assignment problem
  assignment, essential_rows, essential_cols = solve_corresponding_simple_assignment(u,v)
  old_budget_tot = float('inf')
  budget_tot = sum(u) + sum(v)
  while not check_optimal_solution(assignment, essential_rows, essential_cols) and old_budget_tot > budget_tot:
    d = min([u[i]+v[j]-ratings[i][j] for i in range(len(essential_rows)) if not essential_rows[i] for j in range(len(essential_cols)) if not essential_cols[j]])
    min_u = min([u[i] for i in range(len(essential_rows)) if not essential_rows[i]])
    min_v = min([v[j] for j in range(len(essential_cols)) if not essential_cols[j]])
    if min_u > 0:
      m = min(min_u, d)
      for i in range(len(essential_rows)):
        if not essential_rows[i]:  u[i] -= m
      for j in range(len(essential_cols)):
        if essential_cols[j]:  v[j] += m
    else:
      m = min(min_v, d)
      for i in range(len(essential_rows)):
        if essential_rows[i]:  u[i] += m
      for j in range(len(essential_cols)):
        if not essential_cols[j]:  v[j] -= m

    assignment, essential_rows, essential_cols = solve_corresponding_simple_assignment(u,v)

    old_budget_tot = budget_tot
    budget_tot = sum(u) + sum(v)

  return assignment, essential_rows, essential_cols
          
