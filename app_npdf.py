import os
import pandas as pd
import streamlit as st

# Constants
DATA_FILE_PATH = "tool.xlsx"
DATA_FILE_PATH2 = "ddiTable.xlsx"

# Load Data
def load_data(file_path):
    if not os.path.exists(file_path):
        st.error(f"The data file '{file_path}' does not exist.")
        return None
    try:
        return pd.read_excel(file_path)
    except Exception as e:
        st.error(f"An error occurred while loading the data: {e}")
        return None

# Search Drug Interaction
def search_drug_interaction(df, drug1, drug2):
    drug1, drug2 = drug1.strip().lower(), drug2.strip().lower()
    try:
        filtered_data = df[
            (df["Object compound"].str.lower() == drug1)
            & (df["Precipitant compound"].str.lower() == drug2)
        ]
        return filtered_data
    except KeyError as e:
        st.error(f"The expected column '{e}' is missing from the data file.")
        return pd.DataFrame()
    except Exception as e:
        st.error(f"An error occurred during the search: {e}")
        return pd.DataFrame()

# Drug-Drug Interaction Analysis
def ddi_analysis_for_drugs(file_path, drug1, drug2):
    try:
        df = pd.read_excel(file_path)
    except FileNotFoundError:
        return f"Error: The file '{file_path}' was not found.", []
    except Exception as e:
        return f"Error: {str(e)}", []

    try:
        drug1, drug2 = drug1.strip().lower(), drug2.strip().lower()
        df["Compound name"] = df["Compound name"].str.lower()

        drug1_genes = df[df["Compound name"] == drug1]["Interacting gene name"].unique()
        drug2_genes = df[df["Compound name"] == drug2]["Interacting gene name"].unique()

        if not drug1_genes.size:
            return f"Warning: The drug '{drug1}' was not found in the dataset.", []
        if not drug2_genes.size:
            return f"Warning: The drug '{drug2}' was not found in the dataset.", []
    except KeyError as e:
        return f"Error: Column '{e}' not found in the dataframe.", []
    except Exception as e:
        return f"Error: {str(e)}", []

    common_genes = set(drug1_genes) & set(drug2_genes)
    result = []
    victim_perpetrator_results = []
    if common_genes:
        result = [
            f"Common genes found between '{drug1}' and '{drug2}': {', '.join(common_genes)}",
            "Drug-Drug Interaction (DDI) is possible based on these common genes.",
        ]

        interactions = [
            (
                gene,
                (
                    drug1
                    if any(
                        "substrate" in role
                        for role in df[
                            (df["Compound name"] == drug1)
                            & (df["Interacting gene name"] == gene)
                        ]["Compound Gene Relation"].str.lower()
                    )
                    else "-"
                ),
                (
                    drug2
                    if any(
                        "inducer" in role
                        for role in df[
                            (df["Compound name"] == drug2)
                            & (df["Interacting gene name"] == gene)
                        ]["Compound Gene Relation"].str.lower()
                    )
                    else "-"
                ),
                (
                    drug2
                    if any(
                        "inhibitor" in role
                        for role in df[
                            (df["Compound name"] == drug2)
                            & (df["Interacting gene name"] == gene)
                        ]["Compound Gene Relation"].str.lower()
                    )
                    else "-"
                ),
            )
            for gene in common_genes
        ]

        if interactions:
            interaction_df = pd.DataFrame(
                interactions, columns=["Gene", "Substrate", "Inducer", "Inhibitor"]
            )
            result.append(interaction_df)

            for _, row in interaction_df.iterrows():
                if row["Substrate"] == drug1 and (
                    row["Inducer"] == drug2 or row["Inhibitor"] == drug2
                ):
                    victim_perpetrator_results.append(
                        f"'{row['Gene']}', '{drug1}' is likely to be a victim drug."
                    )
                else:
                    victim_perpetrator_results.append(
                        f"'{row['Gene']}', no sufficient information found for '{drug1}' as victim drug and '{drug2}' as perpetrator drug."
                    )

    else:
        result = [
            f"No common genes found between '{drug1}' and '{drug2}'. Cannot report possibility of Drug-Drug Interaction (DDI)."
        ]

    return result, victim_perpetrator_results

# Streamlit UI
st.title("Drug-Drug Interaction Navigator Tool")
st.write("This tool allows you to search for interactions between two drugs. Select drugs to see detailed analysis.")

drug_data = load_data(DATA_FILE_PATH2)
interaction_result = pd.DataFrame()  # Initialize interaction_result as empty DataFrame
ggi_text = []  # Initialize ggi_text as empty list
victim_perpetrator_results = []  # Initialize victim_perpetrator_results as empty list

if drug_data is not None:
    custom_drug_list = sorted(drug_data["Compound name"].dropna().unique().tolist())
    drug1 = st.selectbox("Select the name of the first drug:", custom_drug_list)

    filtered_drug_list = [drug for drug in custom_drug_list if drug != drug1]
    drug2 = st.selectbox("Select the name of the second drug:", filtered_drug_list)

    if st.button("Analyze"):
        df = load_data(DATA_FILE_PATH)
        if df is not None:
            interaction_result = search_drug_interaction(df, drug1, drug2)

            SELECTED_COLUMNS = [
                "Object compound", "Precipitant compound", "Gene Name", "Gene type",
                "Impact on Object concentration", "AUC change (%)", "AUC fold change", 
                "Evidence type", "Reference"
            ]

            SELECTED_COLUMNS_RENAME_MAP = {
                "Object compound": "Victim Drug",
                "Precipitant compound": "Precipitant Drug",
                "Gene Name": "Gene Name",
                "Gene type": "Gene type",
                "Impact on Object concentration": "Impact on Victim Drug concentration",
                "AUC change (%)": "AUC change (%)",
                "AUC fold change": "AUC fold change",
                "Evidence type": "Evidence type",
                "Reference": "Reference"
            }

            st.subheader("DDI Analysis:")
            ggi_text, victim_perpetrator_results = ddi_analysis_for_drugs(DATA_FILE_PATH2, drug1, drug2)
            for line in ggi_text:
                if isinstance(line, str):
                    st.text(line)

            st.subheader("Potential Drug-Drug Interactions Based on their Gene:")
            for line in ggi_text:
                if isinstance(line, pd.DataFrame):
                    st.table(line)

            st.subheader("Victim-Perpetrator Analysis:")
            victim_drug_statement = False
            if victim_perpetrator_results:
                for vp_result in victim_perpetrator_results:
                    if isinstance(vp_result, str):
                        st.text(vp_result)

            st.subheader("Evidence of Documented Drug Interactions:")
            if not interaction_result.empty:
                interaction_result = interaction_result[SELECTED_COLUMNS].rename(columns=SELECTED_COLUMNS_RENAME_MAP)
                st.table(interaction_result)

                # Add victim drug statement below tabulated data only if interactions exist
                for vp_result in victim_perpetrator_results:
                    if "is likely to be a victim drug." in vp_result:
                        victim_drug_statement = True
                if victim_drug_statement:
                    st.markdown(
                        f"\n**Note**: {drug1} is identified as a victim drug and {drug2} is identified as a perpetrator drug in the interactions listed above.")

            else:
                st.write(
                    f"No clinical evidence found between {drug1} and {drug2} to confirm the possibility of DDI."
                )

    if st.button("Reset"):
        st.session_state.clear()
