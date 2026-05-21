import numpy as np
import matplotlib.pyplot as plt

criteria_names = np.array([["flight performance"],
                            ["development risk"],
                            ["cost"],
                            ["logistics"],
                            ["sustainability"]])

options = np.array([["option 1"],
                    ["option 2"],   
                    ["option 3"],
                    ["option 4"]])

weights = np.array([[0.3],
                    [0.2],
                    [0.15],
                    [0.25],
                    [0.1]])

score_table = np.array([[2.0, 2.0, 3.0, 2.0, 2.0],
                         [2.0, 2.0, 2.0, 2.0, 2.0],
                         [2.0, 2.0, 3.0, 3.0, 3.0],
                         [2.0, 1.0, 2.0, 3.0, 2.0]]) # 4 by 5 matrix with the scores for each option and each criterion

Dummy_table = np.array([[0.4, 0.6, 0.4, 0.2, 0.7],
                         [0.7, 0.5, 0.8, 0.8, 0.6],
                         [0.6, 0.5, 0.2, 0.7, 0.1],
                         [1.0, 0.3, 0.9, 0.6, 0.4]])



def calculate_weighted_sum(score_table, weights):
    return score_table @ weights


print(calculate_weighted_sum(score_table, weights))

def plot_sensitivity_analysis(score_table, weights, options):
    weighted_sums = calculate_weighted_sum(score_table, weights)
    plt.bar(options.flatten(), weighted_sums.flatten())
    plt.ylabel('Weighted Sum')
    plt.show()

plot_sensitivity_analysis(score_table, weights, options)
plot_sensitivity_analysis(score_table, np.ones((5,1))*0.2, options)

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

remove_criteria_test(score_table, weights, criteria_names, options)

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

change_weights_test(score_table, weights, criteria_names, options)

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

change_weights_test2(score_table, weights, criteria_names, options)




def change_scores_test(score_table, weights, criteria_names, options):
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
    
    fig_width = max(16, len(x_labels) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, n_options))
    
    for j in range(n_options):
        ax.plot(range(len(x_labels)), perturbations[:, j],
                marker='o', label=options[j][0], color=colors[j], linewidth=2)
    
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Baseline')
    
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, fontsize=9, rotation=45, ha='right', rotation_mode='anchor')
    
    ax.set_ylabel('Weighted Sum')
    ax.legend(loc='lower right')
    ax.set_ylim(0, max(perturbations.flatten()) * 1.1)
    
    plt.tight_layout(pad=2.0)
    plt.show()

change_scores_test(score_table, weights, criteria_names, options)

def change_scores_test(score_table, weights, criteria_names, options):
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
            for delta, sign in [(1.0, '+1.0')]:
                modified_score_table = score_table.copy().astype(float)
                modified_score_table[j, i] += delta
                perturbations.append(calculate_weighted_sum(modified_score_table, weights).flatten())
                x_labels.append(f'{options[j][0]} - {criteria_names[i][0]}\n({sign})')
    
    perturbations = np.array(perturbations)
    
    fig_width = max(16, len(x_labels) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, n_options))
    
    for j in range(n_options):
        ax.plot(range(len(x_labels)), perturbations[:, j],
                marker='o', label=options[j][0], color=colors[j], linewidth=2)
    
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Baseline')
    
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, fontsize=9, rotation=45, ha='right', rotation_mode='anchor')
    
    ax.set_ylabel('Weighted Sum')
    ax.legend(loc='lower right')
    ax.set_ylim(0, max(perturbations.flatten()) * 1.1)
    
    plt.tight_layout(pad=2.0)
    plt.show()

change_scores_test(score_table, weights, criteria_names, options)

def change_scores_test(score_table, weights, criteria_names, options):
    n_criteria = len(criteria_names)
    n_options = len(options)

    baseline = calculate_weighted_sum(score_table, weights).flatten()

    group_labels = []
    plus_perturbations = []
    minus_perturbations = []

    for i in range(n_criteria):
        for j in range(n_options):
            modified_plus = score_table.copy().astype(float)
            modified_plus[j, i] += 1.0
            modified_minus = score_table.copy().astype(float)
            modified_minus[j, i] -= 1.0

            plus_perturbations.append(calculate_weighted_sum(modified_plus, weights).flatten())
            minus_perturbations.append(calculate_weighted_sum(modified_minus, weights).flatten())
            group_labels.append(f'{options[j][0]}\n{criteria_names[i][0]}')

    plus_perturbations = np.array(plus_perturbations)
    minus_perturbations = np.array(minus_perturbations)

    n_groups = len(group_labels)
    x_groups = np.arange(1, n_groups + 1)

    # --- Styling ---
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 10,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': True,
        'grid.color': '#e0e0e0',
        'grid.linestyle': '--',
        'grid.linewidth': 0.7,
        'axes.facecolor': '#fafafa',
        'figure.facecolor': 'white',
    })

    fig_width = max(16, (n_groups + 1) * 0.9)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    ax.set_facecolor('#fafafa')

    palette = plt.cm.tab10(np.linspace(0, 1, n_options))

    for j in range(n_options):
        color = palette[j]

        # Vertical range bars between -1.0 and +1.0 per group
        for g in range(n_groups):
            ax.plot(
                [x_groups[g], x_groups[g]],
                [minus_perturbations[g, j], plus_perturbations[g, j]],
                color=color, linewidth=2.5, alpha=0.3, solid_capstyle='round'
            )

        # Midpoint connecting line
        mid_points = (plus_perturbations[:, j] + minus_perturbations[:, j]) / 2
        x_line = np.concatenate([[0], x_groups])
        y_line = np.concatenate([[baseline[j]], mid_points])
        ax.plot(x_line, y_line, color=color, linewidth=2, label=options[j][0],
                alpha=0.85, zorder=3)

        # +1.0 dots (filled)
        ax.scatter(x_groups, plus_perturbations[:, j],
                   color=color, s=55, zorder=5, edgecolors='white', linewidths=1.2)

        # -1.0 dots (hollow)
        ax.scatter(x_groups, minus_perturbations[:, j],
                   facecolors='white', edgecolors=color, s=55, zorder=5, linewidths=1.8)

        # Baseline dot
        ax.scatter([0], [baseline[j]],
                   color=color, s=70, zorder=5, edgecolors='white', linewidths=1.2)

    # Baseline divider
    ax.axvline(x=0.5, color='#aaaaaa', linestyle=':', linewidth=1.2, alpha=0.8)
    ax.text(0.5, 1.01, 'Baseline', ha='center', va='bottom',
            fontsize=8, color='#888888', transform=ax.get_xaxis_transform())

    ax.set_xticks([0] + list(x_groups))
    ax.set_xticklabels(['Baseline'] + group_labels, fontsize=8.5,
                       rotation=40, ha='right', rotation_mode='anchor')
    ax.set_ylabel('Weighted Sum', fontsize=11, labelpad=10)
    ax.set_ylim(0, 3.5)
    ax.set_xlim(-0.5, n_groups + 0.5)

    # Legend with a subtle frame
    legend = ax.legend(loc='lower right', framealpha=0.9, edgecolor='#dddddd',
                       fontsize=9, title='Options', title_fontsize=9)

    # Add a small legend hint for filled vs hollow
    ax.scatter([], [], color='grey', s=50, edgecolors='white', linewidths=1.2,
               label='+1.0', zorder=5)
    ax.scatter([], [], facecolors='white', edgecolors='grey', s=50, linewidths=1.8,
               label='−1.0', zorder=5)
    ax.legend(loc='lower right', framealpha=0.9, edgecolor='#dddddd',
              fontsize=9, title='Options / perturbation', title_fontsize=9)

    plt.tight_layout(pad=2.5)
    plt.show()

change_scores_test(score_table, weights, criteria_names, options)

# Equal Weights Test