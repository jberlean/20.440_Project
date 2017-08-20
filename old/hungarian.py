# Implement Hungarian algorithm

def solve_simple_assignment(quals, initial = None):
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
    # Determine subset of all transfers using BFS (NOT the DFS proposed in the Kuhn 1955)
    # Each row is only used in at most one transfer
    # When an improving transfer is found, the function returns prematurely

    improving_transfers = []
    nonimproving_transfers = []
    incomplete_transfers = []

    used_rows = [False]*len(assignment)
    
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
          return improving_transfers, nonimproving_transfers
          #continue
        new_j = [new_j for new_j,a in enumerate(assignment[i]) if a][0]
        transfer.append((i,new_j))

        # Search this new column (new_j) for an unassigned row.
        # If none, the transfer is complete (but NOT an improvement)
        # For each unassigned row found, a new transfer to be extended is added to the growing list
        new_rows = [new_i for new_i in range(len(quals)) if quals[new_i][new_j] and not used_rows[new_i]]
        if len(new_rows) == 0:
          nonimproving_transfers.append(transfer)
        else:
          for new_i in new_rows:
            new_incomplete_transfers.append(transfer+[(new_i,new_j)])
            used_rows[new_i] = True
      incomplete_transfers = new_incomplete_transfers

    return (improving_transfers, nonimproving_transfers)

  # Determine initial assignments
  if initial is not None:
    assignment = initial
  else:
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
  def compute_initial(ratings):
    a = [max(row) for row in ratings]
    b = [max([ratings[i][j] for i in range(len(ratings))]) for j in range(len(ratings[0]))]

    if sum(a) <= sum(b):
      return a, [0]*len(b)
    else:
      return [0]*len(a), b

  def solve_corresponding_simple_assignment(ratings, u, v, init_assignment=None):
    # Compute the qualification matrix for the corresponding simple assignment problem,
    # given the budgets assigned to rows (u) and columns (v)
    quals = []
    for i, row in enumerate(ratings):
      q_row = []
      for j, r_ij in enumerate(row):
        q_row.append(u[i]+v[j]==r_ij)
      quals.append(q_row)

    # Solve simple assignment problem corresponding to this qualification matrix
    return solve_simple_assignment(quals, initial=init_assignment)

  def check_optimal_solution(ratings, assignment, u, v):
    assigned_rows = sum([any(row) for row in assignment])
    return assigned_rows==len(ratings)

  # Make ratings a square matrix
  matrix_dim = max(len(ratings), len(ratings[0]))
  add_rows = matrix_dim - len(ratings)
  add_cols = matrix_dim - len(ratings[0])
  ratings_aug = [row + [0]*add_cols for row in ratings] + [[0]*matrix_dim]*add_rows

  # Compute initial budget
  u,v = compute_initial(ratings_aug)

  # Solve simple assignment problem until this solution solves the general assignment problem
  assignment_aug, essential_rows, essential_cols = solve_corresponding_simple_assignment(ratings_aug, u,v)
  old_budget_tot = float('inf')
  budget_tot = sum(u) + sum(v)
  while not check_optimal_solution(ratings_aug, assignment_aug, u, v):
    d = min([u[i]+v[j]-ratings_aug[i][j] for i in range(len(essential_rows)) if not essential_rows[i] for j in range(len(essential_cols)) if not essential_cols[j]])
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

    assignment_aug, essential_rows, essential_cols = solve_corresponding_simple_assignment(ratings_aug, u,v, assignment_aug)

    old_budget_tot = budget_tot
    budget_tot = sum(u) + sum(v)

  # Remove extra rows from assignment matrix
  assignment = [[assignment_aug[i][j] for j in range(len(ratings[0]))] for i in range(len(ratings))]
  #print "Hungarian algorithm complete: sum(weights) = {0}; sum(budget) = {1}".format(sum([ratings[i][j] for i in range(len(assignment)) for j in range(len(assignment[i])) if assignment[i][j]]), sum(u)+sum(v))

  return assignment
          
