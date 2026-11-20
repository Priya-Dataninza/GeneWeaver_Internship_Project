
from textual.app import App, ComposeResult
from textual.widgets import Header, Footer, Static

class GeneWeaverDashboard(App):

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static("GeneWeaver dashboard - results will appear here")
        yield Footer()

def run_dashboard():
    GeneWeaverDashboard().run()

if __name__ == "__main__":
    run_dashboard()