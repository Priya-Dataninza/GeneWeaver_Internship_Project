
from textual.app import App, ComposeResult
from textual.widgets import Header, Footer, Static, DataTable

class GeneWeaverDashboard(App):

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static(
            "GeneWeaver - CRISPR guide matching",
            id="summary"
        )
        yield DataTable(id="results")
        yield Footer()

    def on_mount(self):
        self.query_one(DataTable).add_columns(
            "Record",
            "Start",
            "MatchedDNA",
            "Mismatches",
            "Mismatch positions"
        )

def run_dashboard():
    GeneWeaverDashboard().run()

if __name__ == "__main__":
    run_dashboard()