from django.urls import include, path

from .optimize.OptimizeSequence import SearchMotifView, \
    OptimizeSequenceView

urlpatterns = [
    path('tracks/', include('app.genetic_element.urls')),
    path('species/', include('app.specie.urls')),
    path('contact', include('app.contact.urls')),
    path('workspace/', include('app.workspace.urls')),
    path('constructs/', include('app.construct.urls')),
    path('features/', include('app.parameter.urls')),
    path('optimize_seq/from-sketch', OptimizeSequenceView.as_view(), name='Optimize sequence'),
    path('search-motif', SearchMotifView.as_view(), name='Search motif in sequence')
]
