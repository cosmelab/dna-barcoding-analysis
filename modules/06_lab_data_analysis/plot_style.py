"""
Shared plot styling for lab data analysis.
Consistent colors, white backgrounds, clean layouts.
"""

# Dracula accent colors on white background
COLORS = {
    'purple': '#bd93f9',
    'pink': '#ff79c6',
    'cyan': '#8be9fd',
    'green': '#50fa7b',
    'orange': '#ffb86c',
    'red': '#ff5555',
    'yellow': '#f1fa8c',
    'text': '#2d3436',
    'grid': '#e0e0e0',
}

# Standard layout for all plots
LAYOUT = dict(
    plot_bgcolor='white',
    paper_bgcolor='white',
    font=dict(color=COLORS['text'], size=12),
    title=dict(font=dict(size=18, color=COLORS['text']), x=0.5, xanchor='center'),
    xaxis=dict(
        gridcolor=COLORS['grid'],
        linecolor=COLORS['grid'],
        tickfont=dict(size=11, color=COLORS['text']),
    ),
    yaxis=dict(
        gridcolor=COLORS['grid'],
        linecolor=COLORS['grid'],
        tickfont=dict(size=11, color=COLORS['text']),
    ),
    legend=dict(
        orientation='h',
        yanchor='bottom',
        y=1.02,
        xanchor='center',
        x=0.5,
        bgcolor='rgba(255,255,255,0.8)',
        font=dict(size=11),
    ),
    margin=dict(l=60, r=40, t=80, b=60),
    hovermode='x unified',
    hoverlabel=dict(bgcolor='white', font_size=11),
)


def apply_style(fig):
    """Apply consistent style to a plotly figure."""
    fig.update_layout(**LAYOUT)
    return fig
