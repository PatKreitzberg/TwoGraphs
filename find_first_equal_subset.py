def find_first_equal_subset(X, Y):
    # Map to store {sum: subset_elements} for list X
    sums_x = {0: []}
    for val,label in X:
        new_sums = {}
        for s, subset in sums_x.items():
            new_sum = s + val
            # Only add if we haven't seen this sum before to keep the first subset found
            if new_sum not in sums_x:
                new_sums[new_sum] = subset + [(val, label)]
        sums_x.update(new_sums)

    total_sum = sum([x for x,_ in X]) # Given sum(X) == sum(Y)

    # Map to store {sum: subset_elements} for list Y
    sums_y = {0: []}
    for val,label in Y:
        new_sums = {}
        for s, subset in sums_y.items():
            new_sum = s + val
            if new_sum not in sums_y:
                current_subset_y = subset + [(val,label)]

                # Check for the match immediately
                # 0 < new_sum < total_sum ensures they are proper subsets (A < X, B < Y)
                if 0 < new_sum < total_sum and new_sum in sums_x:
                    return {
                        "subset_A": sums_x[new_sum],
                        "subset_B": current_subset_y,
                        "target_sum": new_sum
                    }

                new_sums[new_sum] = current_subset_y
        sums_y.update(new_sums)

    return None
