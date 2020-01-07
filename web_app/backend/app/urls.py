from django.urls import include, path

from .api.optimize_sequence.OptimizeSequence import OptimizeSequenceFileView, SearchMotifView, \
    OptimizeSequenceSkectherView

urlpatterns = [
    path('tracks', include('app.genetic_element.urls')),
    path('species', include('app.specie.urls')),
    path('contact', include('app.contact.urls')),
    path('workspace/', include('app.workspace.urls')),
    path('constructs', include('app.construct.urls')),
    path('optimize_seq/from-sketch', OptimizeSequenceSkectherView.as_view(), name='Optimize sequence'),
    path('optimize_seq/from-file', OptimizeSequenceFileView.as_view(), name='Optimize sequence'),
    path('search-motif', SearchMotifView.as_view(), name='Search motif in sequence')
]
