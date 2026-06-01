import numpy as np
import matplotlib.pyplot as plt

criteria_names = np.array([["flight performance"],
                            ["development risk"],
                            ["cost"],
                            ["logistics"],
                            ["sustainability"]])

options = np.array([["CWB"],
                    ["CWFC"],   
                    ["FWB"],
                    ["FWFC"]])

options2 = np.array([["Conventional Wing; Batteries (WCB)"],
                    ["Conventional Wing; Fuel Cells (CWFC)"],
                    ["Flying Wing; Batteries (FWB)"],
                    ["Flying Wing; Fuel Cells (FWFC)"]])

weights = np.array([[0.2],
                    [0.2],
                    [0.15],
                    [0.3],
                    [0.15]])

score_table = np.array([[3.0, 2.0, 2.0, 2.0, 2.0],
                        [3.0, 2.0, 2.0, 2.0, 2.0],
                        [2.0, 2.0, 3.0, 3.0, 3.0],
                        [2.0, 1.0, 2.0, 3.0, 2.0]]) # 4 by 5 matrix with the scores for each option and each criterion

Dummy_table = np.array([[0.4, 0.6, 0.4, 0.2, 0.7],
                        [0.7, 0.5, 0.8, 0.8, 0.6],
                        [0.6, 0.5, 0.2, 0.7, 0.1],
                        [1.0, 0.3, 0.9, 0.6, 0.4]])



def calculate_weighted_sum(score_table, weights):
    return score_table @ weights


#print(calculate_weighted_sum(score_table, weights))

def plot_sensitivity_analysis(score_table, weights, options):
    weighted_sums = calculate_weighted_sum(score_table, weights)
    plt.bar(options.flatten(), weighted_sums.flatten())
    plt.ylabel('Weighted Sum')
    plt.ylim(2,3)
    plt.show()

#plot_sensitivity_analysis(score_table, weights, options)
#plot_sensitivity_analysis(score_table, np.ones((5,1))*0.2, options)

def remove_criteria_test(score_table, weights, criteria_names, options):
    n_criteria = len(criteria_names)
    fig, axes = plt.subplots(1, n_criteria, figsize=(4 * n_criteria, 5), sharey=True)
    
    for i, ax in enumerate(axes):
        modified_score_table = np.delete(score_table, i, axis=1)
        modified_weights = np.delete(weights, i, axis=0)
        weighted_sums = calculate_weighted_sum(modified_score_table, modified_weights) / (1 - weights[i])
        ax.bar(options.flatten(), weighted_sums.flatten())
        ax.set_ylabel('Weighted Sum')
        ax.set_title(f'Removed\n{criteria_names[i][0]}', fontsize=10)
        ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.show()

#remove_criteria_test(score_table, weights, criteria_names, options)

def change_weights_test(score_table, weights, criteria_names, options):
    n_criteria = len(criteria_names)
    fig, axes = plt.subplots(2, n_criteria, figsize=(4 * n_criteria, 10), sharey=True)

    for i in range(n_criteria):
        for row, delta in enumerate([-0.05, 0.05]):
            ax = axes[row, i]

            # Adjust weight i by delta, distribute remainder proportionally across others
            modified_weights = weights.copy().astype(float)
            modified_weights[i] += delta

            # Renormalize so weights still sum to 1
            total = modified_weights.sum()
            modified_weights = modified_weights / total

            weighted_sums = calculate_weighted_sum(score_table, modified_weights)

            ax.bar(options.flatten(), weighted_sums.flatten())
            ax.set_ylabel('Weighted Sum')
            sign = '+' if delta > 0 else '-'
            ax.set_title(f'{criteria_names[i][0]}\n(w {sign} 5%)', fontsize=10)
            ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.show()

#change_weights_test(score_table, weights, criteria_names, options)

def change_weights_test2(score_table, weights, criteria_names, options):
    # Normalize weights to ensure they sum to 1
    weights = weights / weights.sum()

    n_criteria = len(criteria_names)
    n_options = len(options)
    
    perturbations = []
    x_labels = []
    
    # Baseline
    baseline = calculate_weighted_sum(score_table, weights).flatten()
    perturbations.append(baseline)
    x_labels.append('Baseline')
    
    # Compute "criteria removed" scores
    removed_scores = []
    for i in range(n_criteria):
        modified_score_table = np.delete(score_table, i, axis=1)
        modified_weights = np.delete(weights, i, axis=0)
        modified_weights = modified_weights / modified_weights.sum()  # Normalize after removal
        weighted_sums = calculate_weighted_sum(modified_score_table, modified_weights).flatten()
        removed_scores.append(weighted_sums)
    removed_scores = np.array(removed_scores)
    
    # For each criterion: -5%, removed, +5%
    removed_x_positions = {}
    for i in range(n_criteria):
        for delta, sign in [(-0.05 * weights[i, 0], '-5%')]:  # Scale delta relative to weight
            modified_weights = weights.copy().astype(float)
            modified_weights[i] += delta
            modified_weights = modified_weights / modified_weights.sum()  # Normalize after perturbation
            perturbations.append(calculate_weighted_sum(score_table, modified_weights).flatten())
            x_labels.append(f'{criteria_names[i][0]}\n({sign})')
        
        removed_x_positions[i] = len(perturbations)
        perturbations.append(removed_scores[i])
        x_labels.append(f'{criteria_names[i][0]}\n(removed)')
        
        for delta, sign in [(0.05 * weights[i, 0], '+5%')]:  # Scale delta relative to weight
            modified_weights = weights.copy().astype(float)
            modified_weights[i] += delta
            modified_weights = modified_weights / modified_weights.sum()  # Normalize after perturbation
            perturbations.append(calculate_weighted_sum(score_table, modified_weights).flatten())
            x_labels.append(f'{criteria_names[i][0]}\n({sign})')
    
    perturbations = np.array(perturbations)
    
    # Scale figure width to number of labels to avoid crowding
    n_labels = len(x_labels)
    fig_width = max(16, n_labels * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, n_options))
    
    for j in range(n_options):
        ax.plot(range(n_labels), perturbations[:, j],
                marker='o', label=options[j][0], color=colors[j], linewidth=2)
    
    # Mark the removed points with X
    for i in range(n_criteria):
        x_pos = removed_x_positions[i]
        for j in range(n_options):
            ax.scatter(x_pos, perturbations[x_pos, j],
                       marker='X', s=150, color=colors[j],
                       edgecolors='black', linewidths=0.8, zorder=5)
    
    # Shade alternating criteria bands
    for i in range(n_criteria):
        if i % 2 == 0:
            ax.axvspan(1 + i * 3 - 0.5, 1 + i * 3 + 2.5, alpha=0.07, color='gray')
    
    ax.scatter([], [], marker='X', s=120, color='gray',
               edgecolors='black', linewidths=0.8, label='Criterion removed')
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Baseline')
    ax.set_xticks(range(n_labels))
    
    # Rotate labels to prevent overlap
    ax.set_xticklabels(x_labels, fontsize=9, rotation=45, ha='right', rotation_mode='anchor')
    
    ax.set_ylabel('Weighted Sum')
    ax.legend(loc='lower right')
    
    # Adjust y-axis limits based on the new scale
    ax.set_ylim(0, 3.5)
    
    # Extra bottom padding so rotated labels aren't clipped
    plt.tight_layout(pad=2.0)
    plt.show()

#change_weights_test2(score_table, weights, criteria_names, options2)

# Equal Weights Test

def change_scores_test(score_table, weights, criteria_names, options, option2):
    n_criteria = len(criteria_names)
    n_options = len(options)
    
    perturbations = []
    x_labels = []
    
    # Baseline
    baseline = calculate_weighted_sum(score_table, weights).flatten()
    perturbations.append(baseline)
    x_labels.append('Baseline')
    
    for i in range(n_criteria):
        for j in range(n_options):
            for delta, sign in [(-1.0, '-1.0'), (1.0, '+1.0')]:
                modified_score_table = score_table.copy().astype(float)
                modified_score_table[j, i] += delta
                perturbations.append(calculate_weighted_sum(modified_score_table, weights).flatten())
                x_labels.append(f'{options[j][0]} - {criteria_names[i][0]}\n({sign})')
    
    perturbations = np.array(perturbations)
    
    fig_width = max(16, len(x_labels) * 1.4)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, n_options))
    
    for j in range(n_options):
        ax.plot(range(len(x_labels)), perturbations[:, j],
                marker='o', label=options2[j][0], color=colors[j], linewidth=2)
    
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Baseline')
    
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, fontsize=7, rotation=45, ha='right', rotation_mode='anchor')
    
    ax.set_ylabel('Weighted Sum')
    ax.legend(loc='lower right')
    ax.set_ylim(0, max(perturbations.flatten()) * 1.1)
    
    plt.tight_layout(pad=2.0, rect=[0.05, 0, 1, 1])
    plt.show()

#change_scores_test(score_table, weights, criteria_names, options, options2)

def change_weights_test2(score_table, weights, criteria_names, options):
    # Normalize weights to ensure they sum to 1
    weights = weights / weights.sum()

    n_criteria = len(criteria_names)
    n_options = len(options)
    
    perturbations = []
    x_labels = []
    
    # Baseline
    baseline = calculate_weighted_sum(score_table, weights).flatten()
    perturbations.append(baseline)
    x_labels.append('Baseline')

    # Equal weights point
    equal_weights = np.ones_like(weights) / len(weights)
    equal_sums = calculate_weighted_sum(score_table, equal_weights).flatten()
    perturbations.append(equal_sums)
    x_labels.append('Equal\nWeights')
    
    # Compute "criteria removed" scores
    removed_scores = []
    for i in range(n_criteria):
        modified_score_table = np.delete(score_table, i, axis=1)
        modified_weights = np.delete(weights, i, axis=0)
        modified_weights = modified_weights / modified_weights.sum()
        weighted_sums = calculate_weighted_sum(modified_score_table, modified_weights).flatten()
        removed_scores.append(weighted_sums)
    removed_scores = np.array(removed_scores)
    
    # For each criterion: -5%, removed, +5%
    removed_x_positions = {}
    for i in range(n_criteria):
        for delta, sign in [(-0.05 * weights[i, 0], '-5%')]:
            modified_weights = weights.copy().astype(float)
            modified_weights[i] += delta
            modified_weights = modified_weights / modified_weights.sum()
            perturbations.append(calculate_weighted_sum(score_table, modified_weights).flatten())
            x_labels.append(f'{criteria_names[i][0]}\n({sign})')
        
        removed_x_positions[i] = len(perturbations)
        perturbations.append(removed_scores[i])
        x_labels.append(f'{criteria_names[i][0]}\n(removed)')
        
        for delta, sign in [(0.05 * weights[i, 0], '+5%')]:
            modified_weights = weights.copy().astype(float)
            modified_weights[i] += delta
            modified_weights = modified_weights / modified_weights.sum()
            perturbations.append(calculate_weighted_sum(score_table, modified_weights).flatten())
            x_labels.append(f'{criteria_names[i][0]}\n({sign})')
    
    perturbations = np.array(perturbations)
    
    n_labels = len(x_labels)
    fig_width = max(16, n_labels * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, n_options))
    
    for j in range(n_options):
        ax.plot(range(n_labels), perturbations[:, j],
                marker='o', label=options[j][0], color=colors[j], linewidth=2)
    
    # Mark the removed points with X
    for i in range(n_criteria):
        x_pos = removed_x_positions[i]
        for j in range(n_options):
            ax.scatter(x_pos, perturbations[x_pos, j],
                       marker='X', s=150, color=colors[j],
                       edgecolors='black', linewidths=0.8, zorder=5)
    
    # Shade alternating criteria bands (offset by 1 to account for equal weights point)
    for i in range(n_criteria):
        if i % 2 == 0:
            ax.axvspan(2 + i * 3 - 0.5, 2 + i * 3 + 2.5, alpha=0.07, color='gray')
    
    ax.scatter([], [], marker='X', s=120, color='gray',
               edgecolors='black', linewidths=0.8, label='Criterion removed')
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Baseline')
    ax.axvline(x=1, color='blue', linestyle=':', linewidth=1, alpha=0.5, label='Equal weights')
    ax.set_xticks(range(n_labels))
    
    ax.set_xticklabels(x_labels, fontsize=11, rotation=45, ha='right', rotation_mode='anchor')
    
    ax.set_ylabel('Weighted Sum', fontsize=14)
    ax.legend(loc='lower left', ncols=2)
    ax.set_ylim(1.6, 2.8)
    
    plt.tight_layout(pad=2.0)
    plt.show()

change_weights_test2(score_table, weights, criteria_names, options2)

def change_scores_test(score_table, weights, criteria_names, options, options2):
    n_criteria = len(criteria_names)
    n_options = len(options)

    perturbations = []
    x_positions = []
    minor_labels = []   # small: option + delta
    major_labels = []   # big: criterion name, placed at group center

    x = 0
    group_centers = []
    group_names = []

    # Baseline
    baseline = calculate_weighted_sum(score_table, weights).flatten()
    perturbations.append(baseline)
    x_positions.append(x)
    minor_labels.append('Baseline')
    x += 1.5  # gap before first group

    for i in range(n_criteria):
        group_start = x
        for j in range(n_options):
            for delta, sign in [(-1.0, '−1'), (1.0, '+1')]:
                modified_score_table = score_table.copy().astype(float)
                modified_score_table[j, i] += delta
                perturbations.append(calculate_weighted_sum(modified_score_table, weights).flatten())
                x_positions.append(x)
                minor_labels.append(f'{options[j][0]}({sign})')
                x += 1
        group_end = x - 1
        group_centers.append((group_start + group_end) / 2)
        group_names.append(criteria_names[i][0])
        x += 1.5  # gap between groups

    perturbations = np.array(perturbations)
    x_positions = np.array(x_positions)

    fig_width = max(16, len(minor_labels) * 1.2)
    fig, ax = plt.subplots(figsize=(fig_width, 7))
    colors = plt.cm.tab10(np.linspace(0, 1, n_options))

    for j in range(n_options):
        ax.plot(x_positions, perturbations[:, j],
                marker='o', label=options2[j][0], color=colors[j], linewidth=2)

    ax.axvline(x=x_positions[0], color='black', linestyle='--', linewidth=1, alpha=0.5, label='Baseline')

    # Minor x-axis labels (option + delta), small font
    ax.set_xticks(x_positions)
    ax.set_xticklabels(minor_labels, fontsize=11, rotation=45, ha='right', rotation_mode='anchor')

    # Major criterion labels above x-axis (below plot area), larger font
    for cx, name in zip(group_centers, group_names):
        ax.text(cx, ax.get_ylim()[1], name,
                ha='center', va='top', fontsize=11, fontweight='bold',
                transform=ax.transData)

    # Shading for alternating criterion groups
    group_edges = []
    x2 = x_positions[1]  # first group starts here
    for i in range(n_criteria):
        left = x2 - 0.5
        right = x2 + n_options * 2 - 0.5
        if i % 2 == 0:
            ax.axvspan(left, right, alpha=0.07, color='gray')
        group_edges.append((left, right))
        x2 += n_options * 2 + 1.5

    ax.set_ylabel('Weighted Sum', fontsize=14)
    ax.legend(loc='lower left', ncols=2)
    #ax.set_ylim(0, max(perturbations.flatten()) * 1.15)
    ax.set_ylim(1.6, 3)

    plt.tight_layout(pad=2.0, rect=[0.05, 0.08, 1, 1])
    plt.show()

change_scores_test(score_table, weights, criteria_names, options, options2)