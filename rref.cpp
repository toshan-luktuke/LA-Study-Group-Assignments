#include <bits/stdc++.h>

using namespace std;

void display(vector<vector<double>> &v);
void display(vector<double> &v);
void display(vector<int> &v);
vector<vector<double>> rref(vector<vector<double>> v);
vector<vector<double>> rref(vector<vector<double>> mat, vector<vector<double>> b);
vector<double> xparticular(vector<vector<double>> &mat);
vector<vector<double>> xnull(vector<vector<double>> &mat);
bool zero_row(vector<double> &v);
bool check_form(vector<vector<double>> &mat);
void exchange_row(vector<vector<double>> &mat, int r1, int r2);
vector<double> leading_zeros(vector<vector<double>> &mat);
void arrange_in_descending(vector<double> &z, vector<vector<double>> &mat);
void arrange_in_descending(vector<double> &z, vector<vector<double>> &mat, vector<vector<double>> &b);
void divide_row(vector<vector<double>> &mat, int val, int row);
void subtract_rows(vector<vector<double>> &mat, int row1, double factor, int row2);
void reduce_col(vector<vector<double>> &mat, int val, int row, int col);
void reduce_col(vector<vector<double>> &mat, vector<vector<double>> &b, int val, int row, int col);
double roundoff(double value, unsigned char prec);

//global vars for rows, cols
int m, n;

int main()
{
    cout << "For a given equation Ax = b" << endl;
    cout << "Input the dimensions of the matrix A (mxn)" << endl;
    cin >> m >> n;
    vector<vector<double>> mat(m, vector<double>(n));
    cout << "Enter the values of the matrix: ";

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> mat[i][j]; //{1,2,3} ->
                              //{4,5,6} :down:
        }
    }
    display(mat);

    cout << "Input the matrix b of Ax = b" << endl;

    vector<vector<double>> b(m, vector<double>(1));
    for (int i = 0; i < m; i++)
    {

        cin >> b[i][0]; //{1,2,3} ->
                        //{4,5,6} :down:
    }

    //get Row-Reduced Echelon Form
    vector<vector<double>> rref_mat = rref(mat);
    vector<vector<double>> aug_mat = rref(mat, b);
    // Get particular solution
    vector<double> particular = xparticular(aug_mat);
    // Get Nullspace
    vector<vector<double>> nullspace = xnull(rref_mat);

    cout << "Row-Reduced Echelon Form: -" << endl;
    display(rref_mat);

    cout << "Particular Solution: -" << endl;
    display(particular);

    cout << "Nullspace: -" << endl;
    display(nullspace);

    cout << "Complete Solution is Particular(Xp) + Nullspace(Xn)" << endl;

    return 0;
}

void display(vector<double> &v)
{
    cout << "[ ";
    for (int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " ";
    }
    cout << "]" << endl;
}

void display(vector<int> &v)
{
    cout << "[ ";
    for (int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " ";
    }
    cout << "]" << endl;
}

void display(vector<vector<double>> &v)
{
    cout << "Matrix is being displayed:-" << endl;
    cout << "----------" << endl;
    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < v[i].size(); j++)
        {

            cout << v[i][j] << " ";
        }
        cout << endl;
    }
    cout << "----------" << endl;
}

vector<double> xparticular(vector<vector<double>> &mat)
{
    //Find all the pivots and their positions
    vector<int> piv_pos;
    vector<double> piv_vals(n, 0);
    for (int i = 0; i < mat.size(); i++)
    {
        for (int j = 0; j < mat[i].size(); j++)
        {
            if (mat[i][j] == 1)
            {
                piv_pos.push_back(j);
                break;
            }
        }
    }

    //Go from bottom to top assigning values
    for (int i = mat.size() - 1; i >= 0; i--)
    {
        // Get each row in reverse-order
        vector<double> co_eff_row = mat[i];
        //Check if this is a zero row
        if (zero_row(co_eff_row))
        {
            continue;
        }
        int preserve = piv_pos.back();
        double w_sum = 0;
        for (int j = 0; j < co_eff_row.size() - 1; j++)
        {
            if (j == preserve)
            {
                continue;
            }
            w_sum = co_eff_row[j] * piv_vals[j];
        }

        piv_vals[piv_pos.back()] = co_eff_row.back() - w_sum;

        piv_pos.pop_back();
    }
    //display(piv_vals);

    return piv_vals;
}

vector<vector<double>> xnull(vector<vector<double>> &mat)
{
    // Identify the pivots and their positions
    // Identify the free variables and their positions
    vector<int> piv_pos;
    vector<int> free_vars(n, 1);
    vector<double> piv_vals(n, 0);
    for (int i = 0; i < mat.size(); i++)
    {
        for (int j = 0; j < mat[i].size(); j++)
        {
            if (mat[i][j] == 1)
            {
                piv_pos.push_back(j);
                break;
            }
        }
    }
    int no_free_vars = n - piv_pos.size();
    //Get the positions of free vars
    for (int i = 0; i < piv_pos.size(); i++)
    {
        free_vars[piv_pos[i]] = 0;
    }

    display(piv_pos);
    display(free_vars);

    //Create a DD vector with dimensions m x no_free_vars
    vector<vector<double>> nullspace(no_free_vars, vector<double>(n));
    int p = 0;
    //Going by each free variable
    //Fill 1's at the nullspace at the index of each corresponding variable
    for (int i = 0; i < free_vars.size(); i++)
    {
        if (free_vars[i] == 0)
        {
            continue;
        }

        // at index of free var it will be 1 and of other free vars it will be 0
        nullspace[p][i] = 1;
        for (int j = 0; j < free_vars.size(); j++)
        {
            if (j == i)
            {
                continue;
            }
            if (free_vars[j] == 1)
            {
                nullspace[p][j] = 0; //Corresponding nullspace is 0
            }
        }
        // Iterate through pivot rows of the main matrice determine the other co-effs
        // Iterate columnwise
        for (int j = 0; j < piv_pos.size(); j++)
        {
            //Kaafi time ke baad ye result aaya
            nullspace[p][piv_pos[j]] = mat[j][i] == 0 ? mat[j][i] : -mat[j][i];
            //The ternary operator is to prevent a -ve 0
        }
        p++;
    }

    // display(nullspace);

    return nullspace;
}

bool zero_row(vector<double> &v)
{
    bool ans = true;
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] != 0)
        {
            ans = false;
            break;
        }
    }
    return ans;
}

vector<vector<double>> rref(vector<vector<double>> mat)
{
    // 1. Format the matrix correctly, all columns with leading zeroes should be at the bottom
    vector<double> no_of_leading_zeros = leading_zeros(mat);
    arrange_in_descending(no_of_leading_zeros, mat);
    // 2. Get the first non-zero element of each row -> convert to 1 -> Subtract from rows below
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (mat[i][j] != 0.0 || mat[i][j] != -0)
            {
                //cout << mat[i][j] << endl;
                // Get the leading value
                double val = mat[i][j];
                // Divide whole row by this value thereby converting leader into 1
                divide_row(mat, val, i);
                // Reduce all the rows above and below this element
                // Current element's column become 0
                reduce_col(mat, val, i, j);
                break;
            }
        }
        // Re-format the matrix
        no_of_leading_zeros = leading_zeros(mat);
        arrange_in_descending(no_of_leading_zeros, mat);
    }

    return mat;
}

vector<vector<double>> rref(vector<vector<double>> mat, vector<vector<double>> b)
{
    // 1. Format the matrix correctly, all columns with leading zeroes should be at the bottom
    vector<double> no_of_leading_zeros = leading_zeros(mat);
    arrange_in_descending(no_of_leading_zeros, mat, b);
    // 2. Get the first non-zero element of each row -> convert to 1 -> Subtract from rows below
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (mat[i][j] != 0.0 || mat[i][j] != -0)
            {
                //cout << mat[i][j] << endl;
                // Get the leading value
                double val = mat[i][j];
                // Divide whole row by this value thereby converting leader into 1
                divide_row(mat, val, i);
                divide_row(b, val, i);
                // Reduce all the rows above and below this element
                // Current element's column become 0
                reduce_col(mat, b, val, i, j);
                break;
            }
        }
        // Re-format the matrix
        no_of_leading_zeros = leading_zeros(mat);
        arrange_in_descending(no_of_leading_zeros, mat, b);
    }

    //return the augmented matrix
    for (int i = 0; i < mat.size(); i++)
    {
        mat[i].push_back(b[i][0]);
    }

    return mat;
}

bool check_form(vector<vector<double>> &mat)
{
    //For every row, the leading element has to be 1
    //Everything above and below the leading element has to be 1
    bool ans = true;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (mat[i][j] != 0 && mat[i][j] == 1)
            {
                //Check above and below if all are 1
                for (int k = 0; k < n; k++)
                {
                    if (mat[i][k] != 0)
                    {
                        ans = false;
                        break;
                    }
                }

                break;
            }
            else if (mat[i][j] != 0 && mat[i][j] != 1)
            {
                ans = false;
                break;
            }
        }
    }

    return ans;
}

void reduce_col(vector<vector<double>> &mat, int val, int row, int col)
{
    for (int i = 0; i < m; i++)
    {
        if (i != row)
        {
            //Subtract all columns by the required values
            double factor = mat[i][col];
            //cout << "Factor: " << factor << endl;
            subtract_rows(mat, i, factor, row);
        }
    }
}

void reduce_col(vector<vector<double>> &mat, vector<vector<double>> &b, int val, int row, int col)
{
    for (int i = 0; i < mat.size(); i++)
    {
        if (i != row)
        {
            //Subtract all columns by the required values
            double factor = mat[i][col];
            //cout << "Factor: " << factor << endl;
            subtract_rows(mat, i, factor, row);
            subtract_rows(b, i, factor, row);
        }
    }
}

void subtract_rows(vector<vector<double>> &mat, int row1, double factor, int row2)
{
    for (int j = 0; j < mat[row1].size(); j++)
    {
        mat[row1][j] -= (factor * mat[row2][j]);
        //cout << "Subtrahend: " << (factor * mat[row2][j]) << "," << mat[row2][j] << " Factor: " << factor << endl;
    }
}

void exchange_row(vector<vector<double>> &mat, int r1, int r2)
{
    vector<double> temp = mat[r2];
    mat[r2] = mat[r1];
    mat[r1] = temp;
}

vector<double> leading_zeros(vector<vector<double>> &mat)
{
    vector<double> no_of_leading_zeros(m, 0);
    //get the no of leading zeros in every row
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (mat[i][j] == 0)
            {
                no_of_leading_zeros[i]++;
            }
            else
            {
                break;
            }
        }
    }
    return no_of_leading_zeros;
}

void arrange_in_descending(vector<double> &z, vector<vector<double>> &mat)
{
    //Selection sort
    int min = z[0];
    for (int i = 0; i < z.size(); i++)
    {
        min = i;
        for (int j = i + 1; j < z.size(); j++)
        {
            if (z[j] < z[min])
            {
                min = j;
            }
        }
        double temp = z[i];
        z[i] = z[min];
        z[min] = temp;
        exchange_row(mat, i, min);
    }

    // for (int i = 0; i < m; i++)
    // {
    //     cout << v[i] << " " << endl;
    // }
}

void arrange_in_descending(vector<double> &z, vector<vector<double>> &mat, vector<vector<double>> &b)
{
    //Selection sort
    int min = z[0];
    for (int i = 0; i < z.size(); i++)
    {
        min = i;
        for (int j = i + 1; j < z.size(); j++)
        {
            if (z[j] < z[min])
            {
                min = j;
            }
        }
        double temp = z[i];
        z[i] = z[min];
        z[min] = temp;
        exchange_row(mat, i, min);
        exchange_row(b, i, min);
    }

    // for (int i = 0; i < m; i++)
    // {
    //     cout << v[i] << " " << endl;
    // }
}

void divide_row(vector<vector<double>> &mat, int val, int row)
{
    for (int j = 0; j < mat[row].size(); j++)
    {
        mat[row][j] = mat[row][j] / val;
        if (mat[row][j] < 0.001 && mat[row][j] > -0.001)
        {
            mat[row][j] = 0.0;
        }
    }
}

double roundoff(double value, unsigned char prec)
{
    double pow_10 = pow(10.0f, (double)prec);
    return round(value * pow_10) / pow_10;
}