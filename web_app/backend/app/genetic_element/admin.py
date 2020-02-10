from django.contrib import admin
from django.contrib.admin import ModelAdmin
from django.db.models import Count
from django.utils.translation import ugettext_lazy as _

from .models import Category, GeneticElement, image_preview, image_list_preview


@admin.register(Category)
class CategoryAdmin(ModelAdmin):
    fieldsets = (
        (_('Information'), {'fields': ('name',)}),
        (_('Important dates'), {'fields': ('created_at', 'updated_at')}),
    )
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('name',),
        }),
    )
    list_display = ('name', 'elements_count', 'created_at')
    search_fields = ('name',)
    ordering = ('-created_at',)

    def get_queryset(self, request):
        return Category.objects.annotate(elements_count=Count('elements'))

    def elements_count(self, obj):
        return obj.elements_count


@admin.register(GeneticElement)
class GeneticElementAdmin(ModelAdmin):
    fieldsets = (
        (_('Image'), {'fields': ('glyph_thumbnail', image_preview)}),
        (_('Information'), {'fields': ('name', 'category', 'role')}),
        (_('Actions'), {'fields': ('visible', 'order', 'default')})
    )
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('name',),
        }),
    )
    date_hierarchy = 'created_at'
    list_display = ('name', 'role', 'order', 'default', image_list_preview, 'category','created_at')
    search_fields = ('name',)
    ordering = ('order',)
    readonly_fields = [image_preview]
