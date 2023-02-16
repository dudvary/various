# Visualizes summary excel sheet (xls sheet, as exported by the Splid App) as Pie Chart and Bar Diagram
# Link to Splid App: https://splid.app/german
#
# Tested for a Splid group of two persons (should work with 2 or more persons)
# with multiple currencies. Script converts currencies to the currency set on the Splid App (here: EUR)
#
# Requirements of exported excel sheet:
# - Splid sheet should contain two sheets named "Zusammenfassung" and "Währungen"
# - Splid needs to be set to German language (to export the excel sheets in German; otherwise the Sheet names will differ)

# Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------
# Modify path to the Splid excel sheet here
filename = "path/to/Zusammenfassung XX.xls"
# ---------------------------------

# FUNCTIONS
def process_splid_files(filename):
    """
    Process Splid Excel sheet (xls format)

    Args:
        filename: str
            filename to Splid generated excel sheet

    Returns:
        df: Pandas DataFrame
            table of expenses (Summary Sheet)
        currencyDict: dict
            currency and its conversion rate to EUR
    """
    # Check sheets of excel file
    xls_file = pd.ExcelFile(filename)
    if 'Währungen' not in xls_file.sheet_names:
        raise NameError('The sheet Währungen does not exist')
    if 'Zusammenfassung' not in xls_file.sheet_names:
        raise NameError('The sheet Zusammenfassung does not exist')
        
    # Load the data from the Splid App summary excel sheet into a pandas DataFrame, skip the first 3 rows
    df = xls_file.parse(header=3,sheet_name='Zusammenfassung') # Load Summary Sheet
    df = df[:-1] # delete last row (contains only summed values)
    dfCurrency = xls_file.parse(header=0,sheet_name='Währungen') # Load Currency Sheet for currency conversion
    
    # Extract conversion rate for currencies from dfCurrency and convert to dictionary -> Currency : ConversionRatio to EUR
    dfCurrency = xls_file.parse(header=0,sheet_name='Währungen')
    # Clean up column
    dfCurrency[dfCurrency.columns[1]] = dfCurrency[dfCurrency.columns[1]].str.replace(r'[()]','') # remove brackets
    dfCurrency[dfCurrency.columns[4]] = dfCurrency[dfCurrency.columns[4]].str.replace(r'[€ ]','') # remove Euro sign and empty space
    dfCurrency[dfCurrency.columns[4]] = dfCurrency[dfCurrency.columns[4]].str.replace(r',','.') # replace , by . for conversion to number
    # Convert conversion rate to number
    dfCurrency[dfCurrency.columns[4]] = pd.to_numeric(dfCurrency[dfCurrency.columns[4]])
    # Convert to dictionary
    currencyDict = dfCurrency.set_index(dfCurrency.columns[1]).to_dict()[dfCurrency.columns[4]]
    
    return df, currencyDict

def pseudonymize_members(df):
    """
    Rename names of members of dataframe (splid Excel sheet)

    Args:
        df: Pandas DataFrame
            table of expenses (Summary Sheet)

    Returns:
        df: Pandas DataFrame
            table of expenses (Summary Sheet) with pseudonymized members names
    """
    member_counter = 1
    for i in range(7,len(df.columns),2):
        new_name = 'member' + str(member_counter)
        df = df.rename(columns={df.columns[i]: new_name})
        member_counter += 1
        
    return df

def get_members(df):
    """
    Determine group members of Splid Excel sheet

    Args:
        df: Pandas DataFrame
            table of expenses (Summary Sheet)

    Returns:
        members: list
            member names
        idx_members_expenses: numpy.ndarray
            boolean of member columns in df
    """
        
    members = list(df.iloc[:, 7::2].columns) # member names
    idx_members = df.columns.isin(members) # their corresponding indices
    members_expenses = list(df.iloc[:, 8::2].columns) # column names of member expenses
    idx_members_expenses = df.columns.isin(members_expenses) # column names of member expenses
    idx_members_val = idx_members_expenses | idx_members  
    
    return members, idx_members_expenses

def clean_df(df,currencyDict, idx_members_val):
    """
    Clean Splid Excel sheet

    Args:
        df: Pandas DataFrame
            table of expenses (Summary Sheet)
        currencyDict: dict
            currency and its conversion rate to EUR
        idx_members_expenses: list
            boolean of member columns in df

    Returns:
        df: Pandas DataFrame
            cleaned table of expenses (Summary Sheet)
    """
    
    # Convert all values to one common currency (here EURO)
    for currency_key, conversion_rate in currencyDict.items():
        idx = (df['Währung']==currency_key)
        if sum(idx)==0:
            print('Currency ' + currency_key + ' was not found in Währungstable!')
        else:
            df.loc[idx,'Währung'] = 'EUR'
            df.loc[idx,'Betrag'] = df.loc[idx,'Betrag']*conversion_rate
            df.loc[idx,idx_members_val] = df.loc[idx,idx_members_val]*conversion_rate

    # Round EURO values to 2 decimals
    df['Betrag'] = round(df['Betrag'],2)
    df.loc[:,idx_members_val] = round(df.loc[:,idx_members_val],2)
    
    # delete payments
    df = delete_payments(df)
    
    return df

def delete_payments(df):
    """
    Delete payment rows form df. 
    Displays number of deleted rows.

    Args:
        df: Pandas DataFrame
            table of expenses (Summary Sheet)

    Returns:
        df: Pandas DataFrame
            table of expenses without payments
    """
    idx_del = (df['Titel']=='Zahlung')
    if sum(idx_del)>0:
        print('Delete ' + str(sum(idx_del)) + ' rows with payments')
    df = df.drop(df[idx_del].index)
    return df

def get_expenses_by_category(df):
    """
    Compute expenses by category.

    Args:
        df: Pandas DataFrame
            table of expenses (Summary Sheet)

    Returns:
        dfCategory: Pandas DataFrame
            table of expenses per Category
    """
    # Only keep relevant columns and change sign of columns -> make positive
    dfRelevant = df[['Kategorie'] + list(df.columns[idx_members_expenses])]
    for i in range(len(members)):
        dfRelevant = dfRelevant.rename(columns={dfRelevant.columns[i+1]: members[i]})
        dfRelevant[members[i]] = dfRelevant[members[i]].abs()
        
    # Create df grouped by categories
    dfCategory = dfRelevant.groupby('Kategorie').sum()
    
    return dfCategory

# VISUALIZATION FUNCTIONS
def plot_pie_chart(df, members):
    """
    plots one pie chart per member.
    Each pie chart shows  expenses per category

    Args:
        df: Pandas DataFrame
            table of expenses (Summary Sheet)
        members: list
            member names        
    Returns:
        None
    """
    # Set Colors
    colors = sns.color_palette("Set3", n_colors=dfCategory.shape[0])
    # Construct graphs
    fig, axs = plt.subplots(1,len(members), figsize=(10,10))
    for m in range(len(members)):
        wedges, texts = axs[m].pie(dfCategory[members[m]],colors=colors,wedgeprops={"edgecolor":"k",'linewidth': 1})
        if m==len(members)-1: # only show legend for first plot
            axs[m].legend(wedges,dfCategory.index,title='Kategorie',loc='center left',bbox_to_anchor=(1,0,0.5,1))
        expenses = dfCategory[members[m]].sum()
        axs[m].set_title(members[m] + '\n' + str(round(expenses,2)) + ' EUR')
    plt.show()
    
def plot_bar_graph(df):
    """
    plots bar graph for each category and member.

    Args:
        df: Pandas DataFrame
            table of expenses (Summary Sheet)     
    Returns:
        None
    """  
    # Set Colors
    colors = sns.color_palette("Set2", n_colors=dfCategory.shape[1])
	# Construct graph
    ax = dfCategory.plot(kind="bar", ylabel='Ausgaben in EUR', color=colors)
    plt.xticks(rotation=45)
    plt.grid(axis="y", which="both", alpha=0.5)
    ax.set_axisbelow(True)
    plt.show()

# Main 
# Process files and load as DataFrame and Dictionary
df, currencyDict = process_splid_files(filename)
# Pseudonymize member names (if desired)
df = pseudonymize_members(df)
# Extract member names and column indices
members, idx_members_expenses = get_members(df)
# Clean DataFrame
df = clean_df(df,currencyDict, idx_members_expenses)
# Get Table grouped by Category
dfCategory = get_expenses_by_category(df)

# Visualize table as pie chart and as bar graph
plot_pie_chart(df, members)
plot_bar_graph(df)