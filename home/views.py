# Create your views here.
from django.shortcuts import render


def home(request):
    context = {
        "project_name": "GeneWeaver",
        "total_sequences": 0,
        "processed_sequences": 0,
        "processing_status": "Ready",
    }

    return render(request, "dashboard/home.html", context)
